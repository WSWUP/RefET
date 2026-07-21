"""Client for the US Bureau of Reclamation AgriMet/Hydromet web service.

Fetches daily and hourly weather data for a station and date range and
returns pandas DataFrames with YEAR / MONTH / DAY (/ HOUR) columns plus one
column per requested pcode, in the native AgriMet units for each pcode.

The service is a CGI script that wraps its output in ``BEGIN DATA`` /
``END DATA`` markers::

    BEGIN DATA
    DATE      ,   FALN MN  ,   FALN MX
    07/01/2015,       66.65,      102.80
    END DATA

Responses are cached on disk when a cache directory is given. Cached
responses for date ranges ending near the present expire after a day by
default, since recent AgriMet data is preliminary and revised; historical
ranges are cached indefinitely. Pass ``max_age`` (seconds) to override.

Data is preliminary and unverified -- see the disclaimer the service returns
with every response.

This module requires pandas (``pip install refet[data]``).
"""
import csv
import datetime as dt
import hashlib
import io
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

DAILY_URL = 'https://www.usbr.gov/pn-bin/daily.pl'
HOURLY_URL = 'https://www.usbr.gov/pn-bin/instant.pl'
# Station metadata: siteid, lat, lon, elevation (m), IANA timezone.
# Linked as "CSV Format" from https://www.usbr.gov/pn/agrimet/location.html
LOCATION_URL = 'https://www.usbr.gov/pn/agrimet/location.csv'
TIMEOUT = 120
RETRIES = 3

# AgriMet publishes a network-wide sensor table rather than per-station
# heights: the RM Young 05103 wind monitor is mounted at 3 m (temperature and
# relative humidity at 2 m).
# https://www.usbr.gov/pn/agrimet/aginfo/sensors.html
DEFAULT_ANEMOMETER_HEIGHT = 3.0

# Cached responses whose range ends within this many days of today are
# treated as current data and expire after CURRENT_MAX_AGE seconds.
CURRENT_WINDOW_DAYS = 30
CURRENT_MAX_AGE = 86400.0

# Gaps are spelled several different ways depending on endpoint and vintage.
MISSING = {'NO RECORD', 'MISSING', 'M', '', '998877', '998877.00', '-9999'}
# The classic Hydromet missing sentinel, sometimes with a fractional part
MISSING_SENTINEL = 998877.0
# Values may carry a trailing quality flag, e.g. "-459.69-" or "57.86^"
QUALITY_FLAG = re.compile(r'[-^]$')

# Plausibility ranges by pcode, in the native AgriMet units for that pcode.
# Readings outside the range are treated as missing (a -459.69 F air
# temperature is a dead sensor, not a measurement). Solar allows a small
# negative for nighttime pyranometer offset. The fallback range only screens
# absurdities for pcodes not listed here.
PLAUSIBLE_RANGES = {
    # Air and dew point temperatures [F]
    'MN': (-60.0, 130.0), 'MX': (-60.0, 130.0), 'MM': (-60.0, 130.0),
    'OB': (-60.0, 130.0), 'YM': (-60.0, 130.0), 'TP': (-60.0, 130.0),
    # Wind speed and gust [mph]
    'UA': (0.0, 150.0), 'WS': (0.0, 150.0), 'WG': (0.0, 150.0),
    # Solar radiation [Langleys (daily) or Langleys/hr (15-min/hourly)]
    'SR': (-10.0, 1500.0), 'SI': (-10.0, 1500.0),
    # Actual vapor pressure [kPa]
    'EA': (0.0, 15.0),
    # Relative humidity [percent]
    'TU': (0.0, 100.0),
}
PLAUSIBLE_FALLBACK = (-100.0, 5000.0)


class AgriMetError(Exception):
    """Raised when the service returns nothing usable."""


def clean_value(cell: str, pcode: str | None = None) -> float | None:
    """One raw cell -> float or None, applying every known missing convention.

    Handles NO RECORD, empty fields, ``m``, the 998877 sentinel, trailing
    quality flags (``-459.69-``, ``57.86^``), and range-checks the result
    against the plausibility range for ``pcode`` (native AgriMet units).
    """
    s = cell.strip()
    if s.upper() in MISSING:
        return None
    s = QUALITY_FLAG.sub('', s).strip()
    if not s or s.upper() in MISSING:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    lo, hi = PLAUSIBLE_RANGES.get((pcode or '').upper(), PLAUSIBLE_FALLBACK)
    if abs(v - MISSING_SENTINEL) < 0.5 or not (lo <= v <= hi):
        return None
    return v


def _url(base, station, start, end, pcodes, hourly=False):
    # Start and end are passed as two successive year/month/day triples; the
    # service reads them positionally, so the order here matters.
    parts = [('station', station)]
    for d in (start, end):
        parts += [('year', d.year), ('month', d.month), ('day', d.day)]
    parts += [('pcode', p) for p in pcodes]
    if hourly:
        parts.append(('print_hourly', 1))
    return base + '?' + urllib.parse.urlencode(parts)


def _cache_fresh(path, max_age):
    if not os.path.exists(path):
        return False
    if max_age is None:
        return True
    return (time.time() - os.path.getmtime(path)) < max_age


def _get(url: str, cache: str | None = None, max_age: float | None = None) -> str:
    """Fetch a URL with retries, reading and writing the disk cache."""
    path = None
    if cache:
        os.makedirs(cache, exist_ok=True)
        key = hashlib.sha1(url.encode()).hexdigest()[:16]
        path = os.path.join(cache, f'agrimet_{key}.txt')
        if _cache_fresh(path, max_age):
            with open(path) as f:
                return f.read()

    last_err = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as r:
                text = r.read().decode('utf-8', errors='replace')
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_err = e
            if attempt < RETRIES - 1:
                time.sleep(2 ** attempt)
    else:
        raise AgriMetError(
            f'request failed after {RETRIES} attempts: {last_err}')

    if path:
        with open(path, 'w') as f:
            f.write(text)
    return text


def _data_block(text):
    """Extract the lines between the BEGIN DATA / END DATA markers."""
    try:
        body = text.split('BEGIN DATA', 1)[1].split('END DATA', 1)[0]
    except IndexError:
        raise AgriMetError(
            'response contained no BEGIN DATA block; the service may be down '
            'or the request may have been rejected')
    return [ln for ln in (line.rstrip() for line in body.splitlines()) if ln.strip()]


def _parse(text, station, pcodes, hourly):
    lines = _data_block(text)
    if len(lines) < 2:
        raise AgriMetError('response contained a header but no data rows')

    # Header cells look like "FALN MN" (daily) or "FALN    OB" (hourly); the
    # pcode is the last whitespace-separated token.
    header = [c.strip() for c in lines[0].split(',')]
    columns = [c.split()[-1] if c.split() else c for c in header[1:]]

    rows, stamps = [], []
    for line in lines[1:]:
        cells = [c.strip() for c in line.split(',')]
        if len(cells) != len(header):
            continue
        stamps.append(cells[0])
        rows.append([clean_value(c, p) for c, p in zip(cells[1:], columns)])

    df = pd.DataFrame(rows, columns=columns, dtype='float64')

    fmt = '%m/%d/%Y %H:%M' if hourly else '%m/%d/%Y'
    when = [dt.datetime.strptime(s, fmt) for s in stamps]
    df.insert(0, 'YEAR', [d.year for d in when])
    df.insert(1, 'MONTH', [d.month for d in when])
    df.insert(2, 'DAY', [d.day for d in when])
    if hourly:
        df.insert(3, 'HOUR', [d.hour for d in when])

    missing_pcodes = [p for p in pcodes if p not in df.columns]
    if missing_pcodes:
        raise AgriMetError(
            f'station {station.upper()} did not return pcodes {missing_pcodes}; '
            f'got {sorted(set(df.columns) - {"YEAR", "MONTH", "DAY", "HOUR"})}')
    if df[list(pcodes)].notna().sum().sum() == 0:
        # An unknown station id is not an HTTP error -- every value just comes
        # back as NO RECORD, so this is where a typo actually surfaces.
        raise AgriMetError(
            f'station {station.upper()} returned no data at all for this range. '
            f'Check the station id and that the period is within its record.')
    return df


def _default_max_age(end):
    """Cache lifetime for a request ending on ``end``: a day for ranges that
    touch the present (preliminary data gets revised), forever otherwise."""
    if end >= dt.date.today() - dt.timedelta(days=CURRENT_WINDOW_DAYS):
        return CURRENT_MAX_AGE
    return None


def _fetch(base, station, start, end, pcodes, cache, hourly, max_age):
    pcodes = list(dict.fromkeys(pcodes))
    if max_age == 'auto':
        max_age = _default_max_age(end)
    # Split on calendar years: one year of hourly data is ~8,800 rows, which
    # the service returns comfortably, and this keeps cache entries reusable.
    frames = []
    year = start.year
    while year <= end.year:
        lo = max(start, dt.date(year, 1, 1))
        hi = min(end, dt.date(year, 12, 31))
        url = _url(base, station, lo, hi, pcodes, hourly=hourly)
        frames.append(_parse(_get(url, cache, max_age), station, pcodes, hourly))
        year += 1
    return pd.concat(frames, ignore_index=True)


def fetch_stations(cache: str | None = None) -> dict:
    """All AgriMet stations from location.csv, keyed by lowercase site id.

    Each value has lat, lon, elev (metres), tz (IANA name) and description.
    Stations with no published elevation are returned with elev=None so the
    caller can ask for it explicitly rather than silently computing ET at
    sea level.
    """
    text = _get(LOCATION_URL, cache, max_age=CURRENT_MAX_AGE if cache else None)
    # The file opens with a stray "Content-disposition:" line before the header
    lines = text.splitlines()
    if lines and lines[0].lower().startswith('content-disposition'):
        lines = lines[1:]
    stations = {}
    for row in csv.DictReader(io.StringIO('\n'.join(lines))):
        site = (row.get('siteid') or '').strip().lower()
        if not site:
            continue

        def num(key):
            try:
                return float((row.get(key) or '').strip())
            except ValueError:
                return None

        # Elevations are metres where a vertical_datum is given; the Great
        # Plains rows carry no elevation at all.
        stations[site] = {
            'description': (row.get('description') or '').strip(),
            'lat': num('latitude'), 'lon': num('longitude'),
            'elev': num('elevation'),
            'elev_units': (row.get('vertical_datum') or '').strip() or 'm',
            'tz': (row.get('timezone') or '').strip(),
            'responsibility': (row.get('responsibility') or '').strip(),
        }
    if not stations:
        raise AgriMetError(f'no stations parsed from {LOCATION_URL}')
    return stations


def station_info(station: str, cache: str | None = None) -> dict:
    """Metadata for one station id, with a helpful error when it is unknown."""
    stations = fetch_stations(cache)
    info = stations.get(station.strip().lower())
    if info is None:
        raise AgriMetError(
            f'unknown AgriMet station "{station}". '
            f'{len(stations)} stations are listed at {LOCATION_URL}')
    if info['responsibility'] == 'great_plains':
        raise AgriMetError(
            f'station "{station}" is a Great Plains site, which the pn-bin '
            f'service does not serve; it needs the gp-bin endpoints instead')
    for key in ('lat', 'lon', 'elev'):
        if info[key] is None:
            raise AgriMetError(
                f'station "{station}" has no published {key}; set it explicitly '
                f'in the config [station] section')
    if info['elev_units'] not in ('m', ''):
        raise AgriMetError(
            f'station "{station}" reports elevation in '
            f'"{info["elev_units"]}", not metres; set elev explicitly')
    return info


def fetch_daily(station, start, end, pcodes, cache=None, max_age='auto'):
    """Daily values for one station.

    Parameters
    ----------
    station : str
        AgriMet station id, e.g. 'FALN'.
    start, end : datetime.date
        Inclusive date range.
    pcodes : list of str
        AgriMet parameter codes. Typical daily pcodes: MN, MX, SR, YM, UA.
    cache : str, optional
        Directory for cached responses.
    max_age : float or None or 'auto', optional
        Cache lifetime in seconds. 'auto' (default) expires ranges that touch
        the last 30 days after a day and keeps older ranges forever.

    Returns
    -------
    pandas.DataFrame
        YEAR, MONTH, DAY columns plus one float column per pcode, in native
        AgriMet units.

    """
    return _fetch(DAILY_URL, station, start, end, pcodes, cache, False, max_age)


def fetch_hourly(station, start, end, pcodes, cache=None, max_age='auto'):
    """Hourly values for one station.

    Parameters are as for :func:`fetch_daily`. Typical hourly pcodes: OB, TP,
    WS, SI. Returns YEAR, MONTH, DAY, HOUR columns plus one column per pcode.
    Hourly rows with missing values are omitted by the service rather than
    flagged, so gaps appear as absent rows.

    """
    return _fetch(HOURLY_URL, station, start, end, pcodes, cache, True, max_age)

"""Minimal client for the US Bureau of Reclamation AgriMet/Hydromet web service.

Fetches daily and hourly weather data for a station and date range and returns
frames shaped the way ``make_report.py`` expects: YEAR / MONTH / DAY (/ HOUR)
columns plus one column per requested pcode.

The service is a CGI script that wraps its output in ``BEGIN DATA`` /
``END DATA`` markers::

    BEGIN DATA
    DATE      ,   FALN MN  ,   FALN MX
    07/01/2015,       66.65,      102.80
    END DATA

Responses are cached on disk when a cache directory is given, so re-running a
build does not re-hit the service.

Data is preliminary and unverified -- see the disclaimer the service returns
with every response.
"""
import csv
import datetime as dt
import hashlib
import io
import os
import re
import urllib.parse
import urllib.request

import pandas as pd

DAILY_URL = 'https://www.usbr.gov/pn-bin/daily.pl'
HOURLY_URL = 'https://www.usbr.gov/pn-bin/instant.pl'
# Station metadata: siteid, lat, lon, elevation (m), IANA timezone.
# Linked as "CSV Format" from https://www.usbr.gov/pn/agrimet/location.html
LOCATION_URL = 'https://www.usbr.gov/pn/agrimet/location.csv'
TIMEOUT = 120

# AgriMet publishes a network-wide sensor table rather than per-station heights:
# the RM Young 05103 wind monitor is mounted at 3 m (temperature/RH at 2 m).
# https://www.usbr.gov/pn/agrimet/aginfo/sensors.html
DEFAULT_ANEMOMETER_HEIGHT = 3.0

# Gaps are spelled several different ways depending on endpoint and vintage.
MISSING = {'NO RECORD', 'MISSING', 'M', '', '998877', '998877.00', '-9999'}
# The classic Hydromet missing sentinel, sometimes with a fractional part
MISSING_SENTINEL = 998877.0
# Values may carry a trailing quality flag, e.g. "-459.69-" or "57.86^"
QUALITY_FLAG = re.compile(r'[-^]$')
# Physically impossible readings that signal a dead sensor (-459.69 F is
# absolute zero; the service emits it rather than a gap marker).
PLAUSIBLE = (-100.0, 5000.0)


def clean_value(cell):
    """One raw cell -> float or None, applying every known missing convention."""
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
    if abs(v - MISSING_SENTINEL) < 0.5 or not (PLAUSIBLE[0] <= v <= PLAUSIBLE[1]):
        return None
    return v


class AgriMetError(Exception):
    """Raised when the service returns nothing usable."""


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


def _get(url, cache=None):
    if cache:
        os.makedirs(cache, exist_ok=True)
        key = hashlib.sha1(url.encode()).hexdigest()[:16]
        path = os.path.join(cache, f'agrimet_{key}.txt')
        if os.path.exists(path):
            with open(path) as f:
                return f.read()
    with urllib.request.urlopen(url, timeout=TIMEOUT) as r:
        text = r.read().decode('utf-8', errors='replace')
    if cache:
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
    return [ln for ln in (l.rstrip() for l in body.splitlines()) if ln.strip()]


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
        rows.append([clean_value(c) for c in cells[1:]])

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


def _fetch(base, station, start, end, pcodes, cache, hourly):
    pcodes = list(dict.fromkeys(pcodes))
    # Split on calendar years: one year of hourly data is ~8,800 rows, which the
    # service returns comfortably, and this keeps cache entries reusable.
    frames = []
    year = start.year
    while year <= end.year:
        lo = max(start, dt.date(year, 1, 1))
        hi = min(end, dt.date(year, 12, 31))
        url = _url(base, station, lo, hi, pcodes, hourly=hourly)
        frames.append(_parse(_get(url, cache), station, pcodes, hourly))
        year += 1
    return pd.concat(frames, ignore_index=True)


def fetch_stations(cache=None):
    """All AgriMet stations from location.csv, keyed by lowercase site id.

    Each value has lat, lon, elev (metres), tz (IANA name) and description.
    Stations with no published elevation are returned with elev=None so the
    caller can ask for it explicitly rather than silently computing ET at
    sea level.
    """
    text = _get(LOCATION_URL, cache)
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
        # Elevations are metres where a vertical_datum is given; the 21
        # Great Plains rows carry no elevation at all.
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


def station_info(station, cache=None):
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


def fetch_daily(station, start, end, pcodes, cache=None):
    """Daily values. Typical pcodes: MN, MX, SR, YM, UA."""
    return _fetch(DAILY_URL, station, start, end, pcodes, cache, hourly=False)


def fetch_hourly(station, start, end, pcodes, cache=None):
    """Hourly values. Typical pcodes: OB, TP, WS, SI."""
    return _fetch(HOURLY_URL, station, start, end, pcodes, cache, hourly=True)


if __name__ == '__main__':
    # Smoke test: python examples/agrimet.py FALN
    import sys
    st = sys.argv[1] if len(sys.argv) > 1 else 'FALN'
    d = fetch_daily(st, dt.date(2015, 7, 1), dt.date(2015, 7, 3),
                    ['MN', 'MX', 'SR', 'YM', 'UA'])
    h = fetch_hourly(st, dt.date(2015, 7, 1), dt.date(2015, 7, 1),
                     ['OB', 'TP', 'WS', 'SI'])
    print(d.to_string(index=False))
    print(h.head().to_string(index=False))

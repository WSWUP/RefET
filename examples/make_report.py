"""Build an interactive HTML report of a year of reference ET for a station.

Reads a TOML station config, computes daily and hourly ASCE reference ET with
this library, optionally validates against Ref-ET .out files, and renders the
result into a single self-contained HTML page.

Usage (from the repository root)::

    pip install -e .              # refet reads its version via importlib.metadata
    pip install pandas pytz
    python examples/make_report.py examples/stations/faln_2015.toml

Every number that appears in the page's prose is a *derived fact* -- see
``derive_facts()``. The template writes them as ``{{placeholder}}`` and this
script substitutes them at build time, so the narrative can never quietly
disagree with the data. An unknown or unused placeholder is a build error, not a
silent miss. The prose around those numbers stays hand-authored in the template.

See examples/README.md for the full contract.
"""
import argparse
import datetime as dt
import json
import os
import re
import sys
import tomllib

import numpy as np
import pandas as pd
import pytz

from refet.daily import Daily
from refet.hourly import Hourly

import agrimet

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']

# Fraction of the seasonal peak hour at which the solar band is called "open".
BAND_THRESHOLD = 0.10
# Standard-time hours counted as night when reporting the after-dark ET share.
NIGHT_HOURS = list(range(20, 24)) + list(range(0, 5))


class ConfigError(Exception):
    """Raised for a station config that cannot produce a report."""


# --------------------------------------------------------------------- config
def load_config(path):
    with open(path, 'rb') as f:
        cfg = tomllib.load(f)
    for key in ('start', 'end'):
        if key not in cfg:
            raise ConfigError(f'{path}: missing required key "{key}"')
    for key in ('start', 'end'):
        if isinstance(cfg[key], dt.datetime):
            cfg[key] = cfg[key].date()
        if not isinstance(cfg[key], dt.date):
            raise ConfigError(
                f'{path}: "{key}" must be a TOML date, e.g. {key} = 2015-01-01')
    if cfg['end'] < cfg['start']:
        raise ConfigError(f'{path}: end ({cfg["end"]}) precedes start ({cfg["start"]})')
    for section in ('daily', 'hourly'):
        if section not in cfg or 'columns' not in cfg[section]:
            raise ConfigError(f'{path}: missing required [{section}.columns]')
    cfg.setdefault('units', {})
    cfg.setdefault('method', 'asce')
    cfg.setdefault('rso_type', None)
    cfg.setdefault('source', {'type': 'csv'})
    cfg['source'].setdefault('type', 'csv')
    if cfg['source']['type'] not in ('csv', 'agrimet'):
        raise ConfigError(f'{path}: [source] type must be "csv" or "agrimet"')
    fill_station_defaults(cfg, path)
    cfg.setdefault('network', '')
    for key in ('title', 'name', 'timezone'):
        if not cfg.get(key):
            raise ConfigError(
                f'{path}: missing required key "{key}" (it is filled in '
                f'automatically only for [source] type = "agrimet")')
    return cfg


def fill_station_defaults(cfg, path):
    """For AgriMet sources, take location/timezone from the published station
    list so a config only has to name a station and a date range.

    Anything set explicitly in the config always wins -- the lookup fills gaps,
    it never overrides.
    """
    cfg.setdefault('station', {})
    if cfg['source']['type'] != 'agrimet':
        return
    station = cfg['source'].get('station')
    if not station:
        raise ConfigError(f'{path}: [source] type="agrimet" requires a station id')
    cache = resolve(cfg['source']['cache']) if cfg['source'].get('cache') else None
    try:
        info = agrimet.station_info(station, cache=cache)
    except agrimet.AgriMetError as e:
        raise ConfigError(str(e))

    for key in ('lat', 'lon', 'elev'):
        cfg['station'].setdefault(key, info[key])
    # AgriMet publishes no per-station anemometer height; the network standard
    # is 3 m. This materially changes ETr, so it is surfaced in the report.
    cfg['station'].setdefault('zw', agrimet.DEFAULT_ANEMOMETER_HEIGHT)
    cfg.setdefault('timezone', info['tz'])
    cfg.setdefault('network', 'AgriMet')
    desc = info['description'] or station.upper()
    # "Fallon, Nevada AgriMet Weather Station" -> "Fallon, Nevada"
    short = re.sub(r'\s*(AgriMet\s*)?Weather Station\s*$', '', desc, flags=re.I) or desc
    cfg.setdefault('title', short)
    cfg.setdefault('name', f'{short} ({station.upper()}) AgriMet')


def make_ticks(dates):
    """Axis tick positions and labels, chosen from the span of the record."""
    n = len(dates)
    parsed = [dt.date.fromisoformat(d) for d in dates]
    multi_year = parsed[0].year != parsed[-1].year
    ticks = []
    if n <= 16:
        step = 1 if n <= 8 else 2
        for i in range(0, n, step):
            ticks.append({'i': i, 'label': parsed[i].strftime('%b %-d')})
    elif n <= 70:
        for i in range(0, n, 7):
            ticks.append({'i': i, 'label': parsed[i].strftime('%b %-d')})
    else:
        # Month starts, thinned so labels never crowd
        starts = [i for i, d in enumerate(parsed) if d.day == 1] or [0]
        every = max(1, round(len(starts) / 12))
        for i in starts[::every]:
            d = parsed[i]
            label = d.strftime('%b')
            if multi_year and (d.month == 1 or i == starts[0]):
                label += d.strftime(' %Y')
            ticks.append({'i': i, 'label': label})
    return ticks


def resolve(path):
    return path if os.path.isabs(path) else os.path.join(ROOT, path)


def station_props(out_path):
    """Scrape station metadata out of a Ref-ET .out file header."""
    props = {}
    with open(out_path) as f:
        lines = f.readlines()
    for line in lines:
        s = line.strip()
        if s.startswith('The anemometer height is'):
            props['zw'] = float(line.split(':')[1].split()[0])
        elif s.startswith('The weather station elevation is'):
            props['elev'] = float(line.split(':')[1].split()[0])
        elif s.startswith('The weather station latitude is'):
            props['lat'] = float(line.split(':')[1].split()[0])
        elif s.startswith('The weather station longitude is'):
            lon = float(line.split(':')[1].split()[0])
            props['lon'] = -lon if 'West' in line else lon
    return props, lines


def read_reference(path):
    """Read the tabular block of a Ref-ET .out file (skipping the units row)."""
    with open(path) as f:
        lines = f.readlines()
    start = [i for i, x in enumerate(lines) if x.startswith(' Mo Day Yr')][0]
    return pd.read_csv(path, sep=r'\s+', index_col=False,
                       skiprows=list(range(start)) + [start + 1])


def station_metadata(cfg):
    """Explicit [station] values win, then the AgriMet station list (already
    merged in by fill_station_defaults), then the daily .out file header."""
    props = {}
    ref = cfg['daily'].get('reference')
    if ref:
        props, _ = station_props(resolve(ref))
    props.update({k: v for k, v in cfg.get('station', {}).items() if v is not None})
    missing = [k for k in ('lat', 'lon', 'elev', 'zw') if k not in props]
    if missing:
        raise ConfigError(
            f'station metadata {missing} not given in [station], not available '
            f'from the AgriMet station list, and no [daily] reference file to '
            f'scrape it from')
    return props


def units_for(cfg, variables):
    """Subset of [units] relevant to one timestep, dropping blanks."""
    return {v: cfg['units'][v] for v in variables
            if cfg['units'].get(v, '').strip()}


def require_columns(df, columns, csv_path):
    missing = {k: v for k, v in columns.items() if v not in df.columns}
    if missing:
        raise ConfigError(
            f'{csv_path}: columns {sorted(missing.values())} not found. '
            f'Available: {sorted(df.columns)}')


# ------------------------------------------------------------------ computing
def load_sources(cfg):
    """Fetch or read the raw daily and hourly frames plus a provenance note."""
    src = cfg['source']
    if src['type'] == 'agrimet':
        station = src.get('station')
        if not station:
            raise ConfigError('[source] type="agrimet" requires a station id')
        cache = resolve(src['cache']) if src.get('cache') else None
        daily = agrimet.fetch_daily(
            station, cfg['start'], cfg['end'],
            [c for c in cfg['daily']['columns'].values()], cache=cache)
        hourly = agrimet.fetch_hourly(
            station, cfg['start'], cfg['end'],
            [c for c in cfg['hourly']['columns'].values()], cache=cache)
        note = (f'AgriMet station <code>{station.upper()}</code>, fetched from '
                f'usbr.gov')
        return daily, hourly, note

    for section in ('daily', 'hourly'):
        if 'csv' not in cfg[section]:
            raise ConfigError(f'[{section}] csv is required for source type "csv"')
    read = lambda p: pd.read_csv(resolve(p), engine='python', na_values='NO RECORD')
    return (read(cfg['daily']['csv']), read(cfg['hourly']['csv']),
            f'local CSV files (<code>{cfg["daily"]["csv"]}</code>)')


def in_range(df, cfg, index_dates):
    """Trim a frame to the configured [start, end] window."""
    lo, hi = cfg['start'].isoformat(), cfg['end'].isoformat()
    keep = [lo <= d <= hi for d in index_dates]
    if not any(keep):
        raise ConfigError(
            f'no records between {lo} and {hi}; the source covers '
            f'{min(index_dates)} to {max(index_dates)}')
    return df[keep]


def build_daily(cfg, props, df):
    c = cfg['daily']
    cols = c['columns']
    require_columns(df, cols, 'daily source')

    df = df.copy()
    df['DATE'] = df[['YEAR', 'MONTH', 'DAY']].apply(
        lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    df = in_range(df, cfg, list(df['DATE']))
    df.set_index('DATE', inplace=True, drop=True)
    df['DOY'] = [int(dt.datetime.strptime(d, '%Y-%m-%d').strftime('%j'))
                 for d in df.index]

    # Raw values go in as-is; refet.units.convert handles the conversion, so the
    # factors live in the library rather than being restated here.
    d = Daily(
        tmin=df[cols['tmin']].values, tmax=df[cols['tmax']].values,
        tdew=df[cols['tdew']].values, rs=df[cols['rs']].values,
        uz=df[cols['uz']].values, zw=props['zw'], elev=props['elev'],
        lat=props['lat'], doy=df['DOY'].values, method=cfg['method'],
        rso_type=cfg['rso_type'],
        input_units=units_for(cfg, ['tmin', 'tmax', 'tdew', 'rs', 'uz']),
    )
    out = pd.DataFrame(index=df.index)
    out['DOY'] = df['DOY']
    out['ETR'], out['ETO'] = d.etr(), d.eto()
    # Ancillary terms come straight off the object, no recomputation
    out['TMIN'], out['TMAX'], out['TMEAN'] = d.tmin, d.tmax, d.tmean
    out['TDEW'], out['RS'], out['RSO'] = d.tdew, d.rs, d.rso
    out['RN'], out['VPD'], out['U2'] = d.rn, d.vpd, d.u2

    if c.get('reference'):
        ref = read_reference(resolve(c['reference']))
        ref['DATE'] = [dt.datetime(y, m, dd).strftime('%Y-%m-%d')
                       for y, m, dd in zip(ref['Yr'], ref['Mo'], ref['Day'])]
        ref.set_index('DATE', inplace=True, drop=True)
        out['ETR_REF'] = ref['ETr']
        out['ETO_REF'] = ref['ETo']
    return out


def build_hourly(cfg, props, df):
    c = cfg['hourly']
    cols = c['columns']
    tz = pytz.timezone(cfg['timezone'])
    require_columns(df, cols, 'hourly source')

    df = df.copy()
    df = in_range(df, cfg, [f'{y:04d}-{m:02d}-{d:02d}' for y, m, d
                            in zip(df['YEAR'], df['MONTH'], df['DAY'])])
    # Raw timestamps are local clock time with DST (localizing drops one hour)
    local = df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: tz.localize(dt.datetime(*x)), axis=1)
    df['DOY'] = local.apply(lambda x: int(x.strftime('%j')))
    df['UTC_HOUR'] = local.apply(lambda x: x.astimezone(pytz.utc).hour)
    # Bin the diurnal axis on this zone's *standard* offset. Clock-local time
    # jumps an hour at the DST boundaries, which smears the solar envelope.
    std_minutes = int(tz.localize(dt.datetime(cfg['start'].year, 1, 1))
                      .utcoffset().total_seconds() // 60)
    std = pytz.FixedOffset(std_minutes)
    cfg['standard_offset'] = f'UTC{std_minutes // 60:+d}'
    df['SOLAR_HOUR'] = local.apply(lambda x: x.astimezone(std).hour)
    df['DATETIME'] = local.apply(
        lambda x: x.astimezone(pytz.utc).strftime('%Y-%m-%d %H:00'))
    df.set_index('DATETIME', inplace=True, drop=True)

    h = Hourly(
        tmean=df[cols['tmean']].values, tdew=df[cols['tdew']].values,
        rs=df[cols['rs']].values, uz=df[cols['uz']].values, zw=props['zw'],
        elev=props['elev'], lat=props['lat'], lon=props['lon'],
        doy=df['DOY'].values, time=df['UTC_HOUR'].values, method=cfg['method'],
        input_units=units_for(cfg, ['tmean', 'tdew', 'rs', 'uz']),
    )
    out = pd.DataFrame(index=df.index)
    out['SOLAR_HOUR'] = df['SOLAR_HOUR']
    # Group the heatmap by *standard-time* calendar date, so its columns line up
    # with the daily series even across a DST shift or a multi-year range.
    out['DATE'] = local.apply(lambda x: x.astimezone(std).strftime('%Y-%m-%d')).values
    out['ETR'], out['ETO'] = h.etr(), h.eto()

    if c.get('reference'):
        ref = read_reference(resolve(c['reference']))
        ref['HOUR'] = (ref['HrMn'] / 100).astype(int)
        ref['DATETIME'] = ref[['Yr', 'Mo', 'Day', 'HOUR']].apply(
            lambda x: tz.localize(dt.datetime(*x))
            .astimezone(pytz.utc).strftime('%Y-%m-%d %H:00'), axis=1)
        ref.set_index('DATETIME', inplace=True, drop=True)
        out['ETR_REF'] = ref['ETr']
        out['ETO_REF'] = ref['ETo']
    return out


def clean(values, nd=3):
    """JSON-safe rounding: NaN and inf become null."""
    return [None if not np.isfinite(v) else round(float(v), nd)
            for v in np.asarray(values, dtype=float)]


def stats(resid):
    resid = resid.dropna()
    return {'n': int(resid.size),
            'max_abs': round(float(resid.abs().max()), 4),
            'rmse': round(float(np.sqrt((resid ** 2).mean())), 4),
            'bias': round(float(resid.mean()), 4)}


# --------------------------------------------------------------------- facts
def solar_band(pivot, columns):
    """First and last standard hour whose mean profile clears the threshold.

    Returns (open, close, width_hours) or (None, None, None) when the window has
    no usable radiation -- a short winter range, or a station reporting no Rs.
    """
    window = pivot[[c for c in columns if c in pivot.columns]]
    if not window.shape[1]:
        return None, None, None
    profile = window.mean(axis=1)
    if not np.isfinite(profile.max()) or profile.max() <= 0:
        return None, None, None
    lit = profile.index[profile >= profile.max() * BAND_THRESHOLD]
    if not len(lit):
        return None, None, None
    lo, hi = int(lit.min()), int(lit.max()) + 1
    return f'{lo:02d}:00', f'{hi:02d}:00', hi - lo


def month_bands(pivot, dates):
    """Solar band per calendar month present in the record.

    Deliberately not hemisphere-aware: "widest" and "narrowest" are measured
    from the data, so a southern-hemisphere station or a three-month range
    describes itself correctly without a special case.
    """
    by_month = {}
    for d in dates:
        by_month.setdefault(d[:7], []).append(d)
    bands = {}
    for key, days in by_month.items():
        # Ignore part-months, whose mean profile is noisy and not comparable
        if len(days) < 20:
            continue
        band = solar_band(pivot, days)
        if band[0]:
            bands[key] = band
    return bands


def derive_facts(cfg, props, daily, hourly, pivot, validation, source_note):
    """Every number the page's prose states, derived from the data.

    Values are pre-formatted strings: the template consumes them verbatim, so
    rounding and wording live here rather than being duplicated in the page.
    """
    etr_total, eto_total = daily['ETR'].sum(), daily['ETO'].sum()
    peak_date = dt.datetime.strptime(daily['ETR'].idxmax(), '%Y-%m-%d')
    northern = props['lat'] >= 0
    start, end = cfg['start'], cfg['end']
    span = (end - start).days + 1

    # Solstice the peak should sit near, by hemisphere. Only meaningful when the
    # record actually spans one -- otherwise the peak is just the range maximum.
    solstice_month = 6 if northern else 12
    solstice = dt.date(peak_date.year, solstice_month, 21)
    offset = (peak_date.date() - solstice).days
    name = 'June' if northern else 'December'
    if not (start <= solstice <= end) or abs(offset) > 21:
        # Too far from the solstice for the comparison to mean anything; a peak
        # 51 days out is a weather event, not a solar one.
        solstice_note = 'the highest single day in the record'
    elif abs(offset) <= 1:
        solstice_note = f'the {name} solstice'
    else:
        solstice_note = (f'{abs(offset)} day{"s" if abs(offset) != 1 else ""} '
                         f'{"after" if offset > 0 else "before"} the {name} solstice')

    bands = month_bands(pivot, list(daily.index))
    if bands:
        widest = max(bands.values(), key=lambda b: b[2])
        narrowest = min(bands.values(), key=lambda b: b[2])
    else:
        # Too short a record to compare months; describe the whole span instead
        whole = solar_band(pivot, list(daily.index))
        widest = narrowest = whole if whole[0] else ('—', '—', 0)

    # Night ET can come out slightly negative at a calm, humid station, which
    # is why this is reported as a share and gated rather than asserted.
    night_frac = (hourly[hourly['SOLAR_HOUR'].isin(NIGHT_HOURS)]['ETR'].sum()
                  / hourly['ETR'].sum())
    monthly = daily.groupby([d[:7] for d in daily.index])['ETR'].sum()
    month_name = lambda key: (MONTHS[int(key[5:7]) - 1]
                              + (f' {key[:4]}' if start.year != end.year else ''))

    if span >= 360:
        period_phrase, period_adjective = 'A full year', 'Annual'
    elif span >= 85:
        period_phrase, period_adjective = f'{span // 30} months', 'Total'
    else:
        period_phrase, period_adjective = f'{span} days', 'Total'

    facts = {
        'title': cfg['title'],
        'station_name': cfg['name'],
        'network': cfg['network'] or 'station',
        'source_note': source_note,
        'start_date': start.strftime('%b %-d, %Y'),
        'end_date': end.strftime('%b %-d, %Y'),
        'period_phrase': period_phrase,
        'period_adjective': period_adjective,
        'standard_offset': cfg['standard_offset'],
        'lat': f'{abs(props["lat"]):.4f}°{"N" if northern else "S"}',
        'lon': f'{abs(props["lon"]):.4f}°{"W" if props["lon"] < 0 else "E"}',
        'elev': f'{props["elev"]:g}',
        'zw': f'{props["zw"]:g}',
        'method': cfg['method'],
        'n_days': f'{int(daily["ETR"].notna().sum()):,}',
        'n_hours': f'{int(hourly["ETR"].notna().sum()):,}',
        'etr_total': f'{etr_total:,.0f}',
        'eto_total': f'{eto_total:,.0f}',
        'etr_eto_ratio': f'{etr_total / eto_total:.2f}',
        'peak_etr': f'{daily["ETR"].max():.1f}',
        'peak_date': peak_date.strftime('%b %-d, %Y'),
        'peak_solstice_note': solstice_note,
        'peak_hourly_etr': f'{hourly["ETR"].max():.2f}',
        'summer_open': widest[0], 'summer_close': widest[1],
        'winter_open': narrowest[0], 'winter_close': narrowest[1],
        'night_share': f'{max(night_frac, 0) * 100:.0f}%',
        'busiest_month': month_name(monthly.idxmax()),
        'quietest_month': month_name(monthly.idxmin()),
    }
    if validation:
        v = validation['daily_etr']
        facts.update({
            'val_rmse': f'{v["rmse"]:.4f}',
            'val_max': f'{v["max_abs"]:.3f}',
            'val_bias': f'{v["bias"]:+.4f}',
            'val_n': f'{v["n"]:,}',
            'val_eto_max': f'{validation["daily_eto"]["max_abs"]:.3f}',
            'val_hourly_rmse': f'{validation["hourly_etr"]["rmse"]:.4f}',
            'val_hourly_max': f'{validation["hourly_etr"]["max_abs"]:.3f}',
            'val_hourly_n': f'{validation["hourly_etr"]["n"]:,}',
        })
    missing = [k for k, v in facts.items() if v is None]
    if missing:
        raise ConfigError(f'facts could not be derived from this data: {missing}')
    return facts, night_frac


# ------------------------------------------------------------------ rendering
PLACEHOLDER = re.compile(r'\{\{\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\}\}')
# Facts the page's own scripts read at runtime, e.g. F_.peak_date
JS_REFERENCE = re.compile(r'\bF_\.([a-zA-Z_][a-zA-Z0-9_]*)\b')
# Optional template region: <!--if validation--> ... <!--endif-->
SECTION = re.compile(r'[ \t]*<!--\s*if\s+(\w+)\s*-->\n?(.*?)'
                     r'[ \t]*<!--\s*endif\s*-->\n?', re.S)


def apply_sections(template, flags):
    """Keep or drop <!--if name--> regions before placeholders are resolved.

    Without this, a report with no reference data would fail the placeholder
    check on validation numbers that legitimately do not exist.
    """
    def sub(m):
        name = m.group(1)
        if name not in flags:
            raise ConfigError(f'template has <!--if {name}--> but no such flag')
        return m.group(2) if flags[name] else ''
    return SECTION.sub(sub, template)


def render(template, facts, data_json, usage_source=None):
    """Substitute {{facts}} and the data payload, strictly.

    Both directions are checked. A placeholder with no fact means the prose
    states something nothing derives; a fact nothing references means the
    narrative silently dropped a number it used to report. Either is a build
    failure -- that is the whole point of the split.
    """
    # Usage is measured against the template *before* optional sections are
    # dropped, so a fact referenced only inside a gated-off region still counts
    # as used rather than looking like dead weight.
    scan = usage_source if usage_source is not None else template
    used = set(PLACEHOLDER.findall(scan)) | set(JS_REFERENCE.findall(scan))
    unknown = sorted(set(PLACEHOLDER.findall(template)) - set(facts))
    if unknown:
        raise ConfigError(
            f'template uses undefined placeholders: {unknown}\n'
            f'Add them to derive_facts() or remove them from the template.')
    unused = sorted(set(facts) - used)
    if unused:
        raise ConfigError(
            f'derived facts never used by the template: {unused}\n'
            f'Reference them as {{{{name}}}} in the prose or F_.name in the page '
            f'script, or drop them from derive_facts().')
    out = PLACEHOLDER.sub(lambda m: facts[m.group(1)], template)
    if '__DATA__' not in out:
        raise ConfigError('template has no __DATA__ placeholder for the payload')
    return out.replace('__DATA__', data_json)


def build_payload(cfg, daily_raw, hourly_raw, source_note):
    props = station_metadata(cfg)
    daily = build_daily(cfg, props, daily_raw)
    hourly = build_hourly(cfg, props, hourly_raw)

    # Diurnal heatmap: mean hourly ETr on a (standard hour x date) grid, with
    # columns reindexed onto the daily series so the two axes stay aligned and
    # a gap in the hourly record shows up as a gap rather than a shift.
    dates = list(daily.index)
    pivot = hourly.pivot_table(index='SOLAR_HOUR', columns='DATE', values='ETR',
                               aggfunc='mean')
    pivot = pivot.reindex(index=range(24), columns=dates)

    validation = None
    if 'ETR_REF' in daily:
        resid = daily['ETR'] - daily['ETR_REF']
        valid = resid.dropna().index
        validation = {
            'daily_etr': stats(resid),
            'daily_eto': stats(daily['ETO'] - daily['ETO_REF']),
            'resid_daily_etr': clean(resid[valid].values, 4),
            'ref_daily_etr': clean(daily.loc[valid, 'ETR_REF'].values),
        }
        if 'ETR_REF' in hourly:
            validation['hourly_etr'] = stats(hourly['ETR'] - hourly['ETR_REF'])

    facts, night_frac = derive_facts(cfg, props, daily, hourly, pivot,
                                     validation, source_note)
    # Flags gate the hand-authored sentences whose claims are data-dependent, so
    # a station that does not behave that way simply does not make the claim.
    flags = {'validation': bool(validation), 'night_et': bool(night_frac >= 0.01)}
    payload = {
        'flags': flags,
        'facts': facts,
        'dates': dates,
        'ticks': make_ticks(dates),
        'daily': {k: clean(daily[k]) for k in
                  ['ETR', 'ETO', 'TMIN', 'TMAX', 'TMEAN', 'TDEW', 'RS', 'RSO',
                   'RN', 'VPD', 'U2']},
        'heatmap': [clean(pivot.loc[hr].values) for hr in range(24)],
        'validation': validation,
    }
    return payload, facts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('config', nargs='?',
                   default=os.path.join(HERE, 'stations', 'faln_2015.toml'),
                   help='TOML station config (default: stations/faln_2015.toml)')
    p.add_argument('--template', default=os.path.join(HERE, 'report_template.html'))
    p.add_argument('--output', help='default: examples/<config name>.html')
    p.add_argument('--json', help='also write the raw payload here')
    args = p.parse_args()

    cfg = load_config(args.config)
    daily_raw, hourly_raw, source_note = load_sources(cfg)
    payload, facts = build_payload(cfg, daily_raw, hourly_raw, source_note)

    output = args.output or os.path.join(
        HERE, os.path.splitext(os.path.basename(args.config))[0] + '_report.html')
    data_json = json.dumps(payload)
    with open(args.template) as f:
        raw_template = f.read()
    template = apply_sections(raw_template, payload['flags'])
    with open(output, 'w') as f:
        f.write(render(template, facts, data_json, usage_source=raw_template))
    if args.json:
        with open(args.json, 'w') as f:
            f.write(data_json)

    print(f"{cfg['name']}, {facts['start_date']} to {facts['end_date']}: "
          f"ETr {facts['etr_total']} mm / ETo {facts['eto_total']} mm over "
          f"{facts['n_days']} days and {facts['n_hours']} hours")
    if payload['validation']:
        for key in ('daily_etr', 'daily_eto', 'hourly_etr'):
            s = payload['validation'].get(key)
            if s:
                print(f"  {key:<11} n={s['n']:<6} rmse={s['rmse']:.4f} "
                      f"max={s['max_abs']:.4f} bias={s['bias']:+.4f} mm")
    else:
        print('  no reference files configured; validation panel omitted')
    print(f'Wrote {output}')


if __name__ == '__main__':
    try:
        main()
    except ConfigError as e:
        sys.exit(f'error: {e}')

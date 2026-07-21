"""Compute the daily and hourly report frames and the page payload.

All reference ET math goes through :class:`refet.Daily` and
:class:`refet.Hourly`; ancillary terms come from their ``results()`` method,
never recomputed here.
"""
import datetime as dt
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

from refet.daily import Daily
from refet.estimate import tdew_from_ea
from refet.hourly import Hourly

from .config import ConfigError, resolve
from .facts import derive_facts


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
        props, _ = station_props(resolve(ref, cfg['_base']))
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


def require_columns(df, columns, source_name):
    missing = {k: v for k, v in columns.items() if v not in df.columns}
    if missing:
        raise ConfigError(
            f'{source_name}: columns {sorted(missing.values())} not found. '
            f'Available: {sorted(df.columns)}')


def vapor_kwargs(df, cols, section):
    """Either a tdew or an ea column feeds the vapor pressure pathway."""
    if 'tdew' in cols:
        return {'tdew': df[cols['tdew']].values}
    if 'ea' in cols:
        return {'ea': df[cols['ea']].values}
    raise ConfigError(f'[{section}.columns] must map either "tdew" or "ea"')


def load_sources(cfg):
    """Fetch or read the raw daily and hourly frames plus a provenance note."""
    src = cfg['source']
    if src['type'] == 'agrimet':
        from refet.io import agrimet

        station = src.get('station')
        if not station:
            raise ConfigError('[source] type="agrimet" requires a station id')
        cache = resolve(src['cache'], cfg['_base']) if src.get('cache') else None
        daily = agrimet.fetch_daily(
            station, cfg['start'], cfg['end'],
            list(cfg['daily']['columns'].values()), cache=cache)
        hourly = agrimet.fetch_hourly(
            station, cfg['start'], cfg['end'],
            list(cfg['hourly']['columns'].values()), cache=cache)
        note = (f'AgriMet station <code>{station.upper()}</code>, fetched from '
                f'usbr.gov')
        return daily, hourly, note

    for section in ('daily', 'hourly'):
        if 'csv' not in cfg[section]:
            raise ConfigError(f'[{section}] csv is required for source type "csv"')

    def read(p):
        return pd.read_csv(resolve(p, cfg['_base']), engine='python',
                           na_values='NO RECORD')

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

    # Raw values go in as-is; refet.units.convert handles the conversion, so
    # the factors live in the library rather than being restated here.
    d = Daily(
        tmin=df[cols['tmin']].values, tmax=df[cols['tmax']].values,
        rs=df[cols['rs']].values, uz=df[cols['uz']].values,
        zw=props['zw'], elev=props['elev'], lat=props['lat'],
        doy=df['DOY'].values, method=cfg['method'], rso_type=cfg['rso_type'],
        input_units=units_for(cfg, ['tmin', 'tmax', 'tdew', 'ea', 'rs', 'uz']),
        **vapor_kwargs(df, cols, 'daily'),
    )
    terms = d.results()
    out = pd.DataFrame(index=df.index)
    out['DOY'] = df['DOY']
    out['ETR'], out['ETO'] = d.etr(), d.eto()
    # Ancillary terms come straight from results(), no recomputation
    out['TMIN'], out['TMAX'], out['TMEAN'] = terms['tmin'], terms['tmax'], terms['tmean']
    out['TDEW'] = terms['tdew'] if 'tdew' in terms else tdew_from_ea(terms['ea'])
    out['RS'], out['RSO'] = terms['rs'], terms['rso']
    out['RN'], out['VPD'], out['U2'] = terms['rn'], terms['vpd'], terms['u2']

    if c.get('reference'):
        ref = read_reference(resolve(c['reference'], cfg['_base']))
        ref['DATE'] = [dt.datetime(y, m, dd).strftime('%Y-%m-%d')
                       for y, m, dd in zip(ref['Yr'], ref['Mo'], ref['Day'])]
        ref.set_index('DATE', inplace=True, drop=True)
        out['ETR_REF'] = ref['ETr']
        out['ETO_REF'] = ref['ETo']
    return out


def build_hourly(cfg, props, df):
    c = cfg['hourly']
    cols = c['columns']
    tz = ZoneInfo(cfg['timezone'])
    require_columns(df, cols, 'hourly source')

    df = df.copy()
    df = in_range(df, cfg, [f'{y:04d}-{m:02d}-{d:02d}' for y, m, d
                            in zip(df['YEAR'], df['MONTH'], df['DAY'])])
    # Raw timestamps are local clock time with DST
    local = df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: dt.datetime(*x, tzinfo=tz), axis=1)
    utc = dt.timezone.utc
    df['DOY'] = local.apply(lambda x: int(x.strftime('%j')))
    df['UTC_HOUR'] = local.apply(lambda x: x.astimezone(utc).hour)
    # Bin the diurnal axis on this zone's *standard* offset. Clock-local time
    # jumps an hour at the DST boundaries, which smears the solar envelope.
    std_minutes = int(dt.datetime(cfg['start'].year, 1, 1, tzinfo=tz)
                      .utcoffset().total_seconds() // 60)
    std = dt.timezone(dt.timedelta(minutes=std_minutes))
    cfg['standard_offset'] = f'UTC{std_minutes // 60:+d}'
    df['SOLAR_HOUR'] = local.apply(lambda x: x.astimezone(std).hour)
    df['DATETIME'] = local.apply(
        lambda x: x.astimezone(utc).strftime('%Y-%m-%d %H:00'))
    df.set_index('DATETIME', inplace=True, drop=True)

    h = Hourly(
        tmean=df[cols['tmean']].values, rs=df[cols['rs']].values,
        uz=df[cols['uz']].values, zw=props['zw'], elev=props['elev'],
        lat=props['lat'], lon=props['lon'], doy=df['DOY'].values,
        time=df['UTC_HOUR'].values, method=cfg['method'],
        input_units=units_for(cfg, ['tmean', 'tdew', 'ea', 'rs', 'uz']),
        **vapor_kwargs(df, cols, 'hourly'),
    )
    out = pd.DataFrame(index=df.index)
    out['SOLAR_HOUR'] = df['SOLAR_HOUR']
    # Group the heatmap by *standard-time* calendar date, so its columns line
    # up with the daily series even across a DST shift or a multi-year range.
    out['DATE'] = local.apply(lambda x: x.astimezone(std).strftime('%Y-%m-%d')).values
    out['ETR'], out['ETO'] = h.etr(), h.eto()

    if c.get('reference'):
        ref = read_reference(resolve(c['reference'], cfg['_base']))
        ref['HOUR'] = (ref['HrMn'] / 100).astype(int)
        ref['DATETIME'] = ref[['Yr', 'Mo', 'Day', 'HOUR']].apply(
            lambda x: dt.datetime(*x, tzinfo=tz)
            .astimezone(utc).strftime('%Y-%m-%d %H:00'), axis=1)
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
    # Flags gate the hand-authored sentences whose claims are data-dependent,
    # so a station that does not behave that way simply does not make the
    # claim.
    flags = {'validation': bool(validation), 'night_et': bool(night_frac >= 0.01)}
    from .render import make_ticks

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

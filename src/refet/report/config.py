"""Station report configuration.

A report is described by a TOML file naming a data source (local CSV files or
an AgriMet station), a date range, the column and unit mappings, and optional
Ref-ET reference output to validate against. Relative paths in a config are
resolved against the directory containing the config file.
"""
import datetime as dt
import os
import re
import tomllib


class ConfigError(Exception):
    """Raised for a station config that cannot produce a report."""


def resolve(path, base):
    """Resolve a config path: absolute paths win, relative ones are taken
    from the directory containing the config file."""
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(base, path))


def load_config(path):
    with open(path, 'rb') as f:
        cfg = tomllib.load(f)
    cfg['_base'] = os.path.dirname(os.path.abspath(path))
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

    Anything set explicitly in the config always wins -- the lookup fills
    gaps, it never overrides.
    """
    cfg.setdefault('station', {})
    if cfg['source']['type'] != 'agrimet':
        return
    from refet.io import agrimet

    station = cfg['source'].get('station')
    if not station:
        raise ConfigError(f'{path}: [source] type="agrimet" requires a station id')
    cache = None
    if cfg['source'].get('cache'):
        cache = resolve(cfg['source']['cache'], cfg['_base'])
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

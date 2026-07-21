"""Derived facts: every number the report's prose states.

Values are pre-formatted strings; the template consumes them verbatim, so
rounding and wording live here rather than being duplicated in the page.
The render step fails the build if the template references a fact that does
not exist, or if a derived fact is never referenced.
"""
import datetime as dt

import numpy as np

from .config import ConfigError

MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']

# Fraction of the seasonal peak hour at which the solar band is called "open".
BAND_THRESHOLD = 0.10
# Standard-time hours counted as night when reporting the after-dark ET share.
NIGHT_HOURS = list(range(20, 24)) + list(range(0, 5))


def month_day(d):
    """'Jul 1' without strftime's platform-dependent day modifier.

    ``%-d`` is glibc-only and raises on native Windows Python, so day numbers
    are formatted manually.
    """
    return f'{d.strftime("%b")} {d.day}'


def month_day_year(d):
    """'Jul 1, 2015', Windows-safe (see month_day)."""
    return f'{d.strftime("%b")} {d.day}, {d.year}'


def solar_band(pivot, columns):
    """First and last standard hour whose mean profile clears the threshold.

    Returns (open, close, width_hours) or (None, None, None) when the window
    has no usable radiation -- a short winter range, or a station reporting
    no Rs.
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
    """Every number the page's prose states, derived from the data."""
    etr_total, eto_total = daily['ETR'].sum(), daily['ETO'].sum()
    peak_date = dt.datetime.strptime(daily['ETR'].idxmax(), '%Y-%m-%d')
    northern = props['lat'] >= 0
    start, end = cfg['start'], cfg['end']
    span = (end - start).days + 1

    # Solstice the peak should sit near, by hemisphere. Only meaningful when
    # the record actually spans one -- otherwise the peak is just the range
    # maximum.
    solstice_month = 6 if northern else 12
    solstice = dt.date(peak_date.year, solstice_month, 21)
    offset = (peak_date.date() - solstice).days
    name = 'June' if northern else 'December'
    if not (start <= solstice <= end) or abs(offset) > 21:
        # Too far from the solstice for the comparison to mean anything; a
        # peak 51 days out is a weather event, not a solar one.
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
        # Too short a record to compare months; describe the whole span
        whole = solar_band(pivot, list(daily.index))
        widest = narrowest = whole if whole[0] else ('n/a', 'n/a', 0)

    # Night ET can come out slightly negative at a calm, humid station, which
    # is why this is reported as a share and gated rather than asserted.
    night_frac = (hourly[hourly['SOLAR_HOUR'].isin(NIGHT_HOURS)]['ETR'].sum()
                  / hourly['ETR'].sum())
    monthly = daily.groupby([d[:7] for d in daily.index])['ETR'].sum()

    def month_name(key):
        return (MONTHS[int(key[5:7]) - 1]
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
        'start_date': month_day_year(start),
        'end_date': month_day_year(end),
        'period_phrase': period_phrase,
        'period_adjective': period_adjective,
        'standard_offset': cfg['standard_offset'],
        'lat': f'{abs(props["lat"]):.4f}\N{DEGREE SIGN}{"N" if northern else "S"}',
        'lon': f'{abs(props["lon"]):.4f}\N{DEGREE SIGN}{"W" if props["lon"] < 0 else "E"}',
        'elev': f'{props["elev"]:g}',
        'zw': f'{props["zw"]:g}',
        'method': cfg['method'],
        'n_days': f'{int(daily["ETR"].notna().sum()):,}',
        'n_hours': f'{int(hourly["ETR"].notna().sum()):,}',
        'etr_total': f'{etr_total:,.0f}',
        'eto_total': f'{eto_total:,.0f}',
        'etr_eto_ratio': f'{etr_total / eto_total:.2f}',
        'peak_etr': f'{daily["ETR"].max():.1f}',
        'peak_date': month_day_year(peak_date),
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

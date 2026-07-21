"""Regenerate the FALN 2015 test CSVs from the AgriMet service.

Uses the packaged client (refet.io.agrimet), which handles the BEGIN/END DATA
framing, missing-value conventions, and per-pcode plausibility checks. Run
from this directory::

    python setup_test_data.py

Note the output is written by pandas rather than copied verbatim, so gap
markers appear as empty fields instead of the service's NO RECORD strings.
The test suite reads both.
"""
import datetime as dt

from refet.io import agrimet

station = 'FALN'
year = 2015

print('Retrieving daily data')
daily = agrimet.fetch_daily(
    station, dt.date(year, 1, 1), dt.date(year, 12, 31),
    ['MN', 'MX', 'SR', 'YM', 'UA', 'ETRS', 'ETOS'])
daily_csv = f'{station}_Agrimet_daily_raw_{year}.csv'
daily.to_csv(daily_csv, index=False)
print(f'Wrote {daily_csv} ({len(daily)} rows)')

print('Retrieving hourly data')
hourly = agrimet.fetch_hourly(
    station, dt.date(year, 1, 1), dt.date(year, 12, 31),
    ['OB', 'TP', 'WS', 'SI'])
hourly_csv = f'{station}_Agrimet_hourly_raw_{year}.csv'
hourly.to_csv(hourly_csv, index=False)
print(f'Wrote {hourly_csv} ({len(hourly)} rows)')

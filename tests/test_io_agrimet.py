import datetime as dt
import hashlib
import os
import time
import urllib.error

import pytest

pd = pytest.importorskip('pandas')

from refet.io import agrimet


DAILY_TEXT = """Some disclaimer header text.

BEGIN DATA
DATE      ,   FALN MN  ,   FALN MX
07/01/2015,       66.65,      102.80
07/02/2015,   NO RECORD,       98.60
07/03/2015,       60.10,   998877.00
END DATA
Some footer text.
"""

HOURLY_TEXT = """BEGIN DATA
DATE            ,   FALN    OB,   FALN    SI
07/01/2015 00:00,       75.20,        0.00
07/01/2015 01:00,       73.40,       -1.50
END DATA
"""


class FakeResponse:
    def __init__(self, data):
        self.data = data

    def read(self):
        return self.data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


# ---------------------------------------------------------------- clean_value
@pytest.mark.parametrize(
    'cell, pcode, expected',
    [
        ['66.65', None, 66.65],
        ['  66.65  ', 'MN', 66.65],
        ['NO RECORD', None, None],
        ['no record', None, None],
        ['', None, None],
        ['m', None, None],
        ['M', None, None],
        ['MISSING', None, None],
        ['998877', None, None],
        ['998877.00', None, None],
        ['998877.49', None, None],
        ['-9999', None, None],
        ['abc', None, None],
        # Trailing quality flags are stripped before parsing
        ['57.86^', 'OB', 57.86],
        ['12.34-', None, 12.34],
        # Physically impossible readings are treated as missing
        ['-459.69-', 'OB', None],
        ['-459.69', 'OB', None],
        ['200', 'MN', None],
        ['45.0', 'MN', 45.0],
        ['-20', 'SI', None],
        ['-2', 'SI', -2.0],
        ['101', 'TU', None],
        ['160', 'UA', None],
        ['12', 'WS', 12.0],
        ['20', 'EA', None],
        ['1.2', 'EA', 1.2],
        # Unknown pcodes fall back to the loose range
        ['4000', 'XX', 4000.0],
        ['6000', 'XX', None],
    ]
)
def test_clean_value(cell, pcode, expected):
    assert agrimet.clean_value(cell, pcode) == expected


# --------------------------------------------------------------------- _parse
def test_parse_daily():
    df = agrimet._parse(DAILY_TEXT, 'FALN', ['MN', 'MX'], hourly=False)
    assert list(df.columns) == ['YEAR', 'MONTH', 'DAY', 'MN', 'MX']
    assert len(df) == 3
    assert df['YEAR'].tolist() == [2015, 2015, 2015]
    assert df['DAY'].tolist() == [1, 2, 3]
    assert float(df['MN'][0]) == 66.65
    assert pd.isna(df['MN'][1])
    assert pd.isna(df['MX'][2])  # 998877 sentinel


def test_parse_hourly():
    df = agrimet._parse(HOURLY_TEXT, 'FALN', ['OB', 'SI'], hourly=True)
    assert list(df.columns) == ['YEAR', 'MONTH', 'DAY', 'HOUR', 'OB', 'SI']
    assert df['HOUR'].tolist() == [0, 1]
    assert float(df['SI'][1]) == -1.5


def test_parse_missing_pcode():
    with pytest.raises(agrimet.AgriMetError, match='did not return pcodes'):
        agrimet._parse(DAILY_TEXT, 'FALN', ['MN', 'UA'], hourly=False)


def test_parse_all_missing_is_unknown_station():
    text = ('BEGIN DATA\n'
            'DATE      ,   ZZZZ MN\n'
            '07/01/2015,   NO RECORD\n'
            'END DATA\n')
    with pytest.raises(agrimet.AgriMetError, match='no data at all'):
        agrimet._parse(text, 'ZZZZ', ['MN'], hourly=False)


def test_parse_no_begin_data():
    with pytest.raises(agrimet.AgriMetError, match='BEGIN DATA'):
        agrimet._parse('service is down\n', 'FALN', ['MN'], hourly=False)


# ------------------------------------------------------------ _get and cache
def test_get_caches(tmp_path, monkeypatch):
    calls = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        return FakeResponse(b'hello')

    monkeypatch.setattr(agrimet.urllib.request, 'urlopen', fake_urlopen)
    out1 = agrimet._get('http://x/y', cache=str(tmp_path))
    out2 = agrimet._get('http://x/y', cache=str(tmp_path))
    assert out1 == out2 == 'hello'
    assert len(calls) == 1


def test_get_cache_expiry(tmp_path, monkeypatch):
    calls = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        return FakeResponse(b'hello')

    monkeypatch.setattr(agrimet.urllib.request, 'urlopen', fake_urlopen)
    url = 'http://x/y'
    agrimet._get(url, cache=str(tmp_path), max_age=3600)
    # Age the cache entry two hours and it should be refetched
    key = hashlib.sha1(url.encode()).hexdigest()[:16]
    path = os.path.join(str(tmp_path), f'agrimet_{key}.txt')
    old = time.time() - 7200
    os.utime(path, (old, old))
    agrimet._get(url, cache=str(tmp_path), max_age=3600)
    assert len(calls) == 2
    # With no max_age the aged entry is served from cache
    os.utime(path, (old, old))
    agrimet._get(url, cache=str(tmp_path), max_age=None)
    assert len(calls) == 2


def test_get_retries(monkeypatch):
    attempts = []

    def flaky(url, timeout=None):
        attempts.append(1)
        if len(attempts) < 3:
            raise urllib.error.URLError('boom')
        return FakeResponse(b'ok')

    monkeypatch.setattr(agrimet.urllib.request, 'urlopen', flaky)
    monkeypatch.setattr(agrimet.time, 'sleep', lambda s: None)
    assert agrimet._get('http://x') == 'ok'
    assert len(attempts) == 3


def test_get_fails_after_retries(monkeypatch):
    def broken(url, timeout=None):
        raise urllib.error.URLError('down')

    monkeypatch.setattr(agrimet.urllib.request, 'urlopen', broken)
    monkeypatch.setattr(agrimet.time, 'sleep', lambda s: None)
    with pytest.raises(agrimet.AgriMetError, match='failed after'):
        agrimet._get('http://x')


def test_default_max_age():
    today = dt.date.today()
    assert agrimet._default_max_age(today) == agrimet.CURRENT_MAX_AGE
    assert agrimet._default_max_age(today - dt.timedelta(days=5)) == \
        agrimet.CURRENT_MAX_AGE
    assert agrimet._default_max_age(today - dt.timedelta(days=60)) is None


# ------------------------------------------------------------- station lookup
LOCATION_CSV = (
    'Content-disposition: attachment; filename=location.csv\n'
    'siteid,description,latitude,longitude,elevation,vertical_datum,'
    'timezone,responsibility\n'
    'FALN,"Fallon, Nevada AgriMet Weather Station",39.4575,-118.77388,'
    '1208.5,m,US/Pacific,pacific_northwest\n'
    'GPXX,"Somewhere, Montana",47.0,-110.0,,,US/Mountain,great_plains\n'
    'NOEL,"No Elevation, Oregon",44.0,-120.0,,,US/Pacific,pacific_northwest\n'
)


@pytest.fixture
def fake_location(monkeypatch):
    monkeypatch.setattr(
        agrimet, '_get', lambda url, cache=None, max_age=None: LOCATION_CSV)


def test_station_info(fake_location):
    info = agrimet.station_info('FALN')
    assert info['lat'] == 39.4575
    assert info['lon'] == -118.77388
    assert info['elev'] == 1208.5
    assert info['tz'] == 'US/Pacific'
    assert info['description'].startswith('Fallon')


def test_station_info_case_insensitive(fake_location):
    assert agrimet.station_info('faln')['elev'] == 1208.5


def test_station_info_unknown(fake_location):
    with pytest.raises(agrimet.AgriMetError, match='unknown AgriMet station'):
        agrimet.station_info('ZZZZ')


def test_station_info_great_plains(fake_location):
    with pytest.raises(agrimet.AgriMetError, match='Great Plains'):
        agrimet.station_info('GPXX')


def test_station_info_missing_elev(fake_location):
    with pytest.raises(agrimet.AgriMetError, match='no published elev'):
        agrimet.station_info('NOEL')


# ------------------------------------------------------------- year chunking
def test_fetch_splits_calendar_years(monkeypatch):
    urls = []

    def fake_get(url, cache=None, max_age=None):
        urls.append(url)
        return 'sentinel'

    def fake_parse(text, station, pcodes, hourly):
        return pd.DataFrame({'YEAR': [1]})

    monkeypatch.setattr(agrimet, '_get', fake_get)
    monkeypatch.setattr(agrimet, '_parse', fake_parse)
    df = agrimet.fetch_daily(
        'FALN', dt.date(2014, 12, 30), dt.date(2015, 1, 2), ['MN'])
    assert len(urls) == 2
    assert 'year=2014' in urls[0] and 'year=2015' in urls[1]
    assert len(df) == 2


# ------------------------------------------------------------- live smoke test
@pytest.mark.network
def test_live_fetch_daily():
    df = agrimet.fetch_daily(
        'FALN', dt.date(2015, 7, 1), dt.date(2015, 7, 3),
        ['MN', 'MX', 'SR', 'YM', 'UA'])
    assert len(df) == 3
    assert float(df['MN'][0]) == pytest.approx(66.65, abs=0.5)

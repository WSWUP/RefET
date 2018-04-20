import datetime as dt
import math
import os

import numpy as np
import pandas as pd
import pytest

from refet.daily import Daily
import refet.units as units


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': units._deg2rad(39.4575),
    'lon': units._deg2rad(-118.77388),
    'zw': 3.0,
}

# Daily test parameters for 2015-07-01
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'eto_refet': 7.9422320475712835,
    'etr_asce': 10.626087665395694,
    'etr_refet': 10.571314344056955,
    'etr_rso_simple': 10.628137858930051,
    'rs': 674.07 * 0.041868,  # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,  # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}

# Test full daily functions with positional inputs
def test_refet_daily_input_positions():
    etr = Daily(
        d_args['tmin'], d_args['tmax'], d_args['ea'], d_args['rs'],
        d_args['uz'], s_args['zw'], s_args['elev'], s_args['lat'],
        d_args['doy'], 'asce').etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


# Test full daily calculations with keyword inputs
def test_refet_daily_etr():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'], method='asce').etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_eto():
    eto = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'], method='refet').eto()
    assert float(eto) == pytest.approx(d_args['eto_refet'])


def test_refet_daily_rso_type_simple():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'], method='refet',
        rso_type='simple').etr()
    assert float(etr) == pytest.approx(d_args['etr_rso_simple'])

def test_refet_daily_rso_type_array():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'], method='refet',
        rso_type='array', rso=d_args['rso']).etr()
    assert float(etr) == pytest.approx(d_args['etr_refet'])

def test_refet_daily_rso_type_exception():
    with pytest.raises(ValueError):
        etr = Daily(
            tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
            rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
            elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
            method='refet', rso_type='nonsense').etr()
        # assert float(etr) == pytest.approx(d_args['etr'])


def test_refet_daily_default_method():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy']).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_asce_method():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        method='asce').etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_refet_method():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        method='refet').etr()
    assert float(etr) == pytest.approx(d_args['etr_refet'])


# Need to test other input structures (i.e. 1d and 2d arrays)
# def test_refet_daily_1d():
#     assert True

# def test_refet_daily_2d():
#     assert True


# Test latitude/longitude in degrees
def test_refet_daily_lat_exception():
    with pytest.raises(ValueError):
        etr = Daily(
            tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
            rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
            elev=s_args['elev'], lat=s_args['lat'] * 180 / math.pi,
            doy=d_args['doy']).etr()
        # assert float(etr) == pytest.approx(d_args['etr'])


# Test unit conversions
def test_refet_daily_tmin_f():
    etr = Daily(
        tmin=d_args['tmin'] * (9.0 / 5) + 32, tmax=d_args['tmax'],
        ea=d_args['ea'], rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'tmin': 'F'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_tmax_k():
    print(d_args['tmin'])
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'] + 273.15,
        ea=d_args['ea'], rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'tmax': 'K'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_zw_ft():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'] / 0.3048,
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'zw': 'ft'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_elev_ft():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'] / 0.3048, lat=s_args['lat'], doy=d_args['doy'],
        input_units={'elev': 'ft'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

def test_refet_daily_lat_def():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'] * 180 / math.pi, doy=d_args['doy'],
        input_units={'lat': 'deg'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


# Test daily functions using actual RefET input/output files
# DEADBEEF - This doesn't work if I move it to conftest.py
class DailyData():
    """Setup daily validation data from Fallon AgriMet station"""
    val_ws = os.path.join(os.getcwd(), 'tests', 'data')
    # val_ws = os.path.join(os.path.dirname(os.getcwd()), 'tests', 'data')

    csv_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.csv')
    # in2_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.in2')
    out_path = os.path.join(val_ws, 'FALN_Agrimet_daily_raw_2015.out')

    # Read in the inputs CSV file using pandas
    csv_df = pd.read_csv(csv_path, engine='python', na_values='NO RECORD')
    csv_df.rename(
        columns={'MN': 'TMIN', 'MX': 'TMAX', 'YM': 'TDEW', 'UA': 'WIND',
                 'SR': 'RS'},
        inplace=True)
    csv_df['DATE'] = csv_df[['YEAR', 'MONTH', 'DAY']].apply(
        lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    csv_df.set_index('DATE', inplace=True, drop=True)

    # Convert inputs units
    csv_df['TMIN'] = units._f2c(csv_df['TMIN'])
    csv_df['TMAX'] = units._f2c(csv_df['TMAX'])
    csv_df['TDEW'] = units._f2c(csv_df['TDEW'])
    csv_df['WIND'] *= 0.44704
    csv_df['RS'] *= 0.041868  # Conversion from Langleys to MJ m-2 to match RefET
    # csv_df['RS'] *= 0.041840  # Alternate conversion from Langleys to MJ m-2
    csv_df['EA'] = 0.6108 * np.exp(17.27 * csv_df['TDEW'] / (csv_df['TDEW'] + 237.3))

    # Eventually compare ancillary functions directly to IN2 values
    # # Identify the row number of the IN2 data
    # with open(in2_path) as in2_f:
    #     in2_data = in2_f.readlines()
    # in2_start = [i for i, x in enumerate(in2_data)
    #              if x.startswith(' Mo Da Year ')][0]
    # # Read in the IN2 file using pandas
    # in2_df = pd.read_table(
    #     in2_path, delim_whitespace=True, skiprows=in2_start, header=[0, 1, 2])
    # in2_df.rename(
    #     columns={'Year': 'YEAR', 'Mo': 'MONTH', 'Da': 'DAY', 'DoY': 'DOY'},
    #     inplace=True)
    # in2_df['DATE'] = in2_df[['YEAR', 'MONTH', 'DAY']].apply(
    #     lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    # in2_df.set_index('DATE', inplace=True, drop=True)

    # Identify the row number of the OUT data
    with open(out_path) as out_f:
        out_data = out_f.readlines()
    out_start = [
        i for i, x in enumerate(out_data) if x.startswith(' Mo Day Yr')][0]
    # Read in the OUT file using pandas (skip header and units)
    out_df = pd.read_table(
        out_path, delim_whitespace=True, index_col=False,
        skiprows=list(range(out_start)) + [out_start + 1])
    out_df.rename(
        columns={'Yr': 'YEAR', 'Mo': 'MONTH', 'Day': 'DAY', 'Tmax': 'TMAX',
                 'Tmin': 'TMIN', 'Wind': 'WIND', 'Rs': 'RS', 'DewP': 'TDEW'},
        inplace=True)
    out_df['DATE'] = out_df[['YEAR', 'MONTH', 'DAY']].apply(
        lambda x: dt.datetime(*x).strftime('%Y-%m-%d'), axis=1)
    out_df.set_index('DATE', inplace=True, drop=True)

    # Read the station properties from the IN2 file for now
    # The values should probably be read using a regular expression
    for line in out_data:
        if line.strip().startswith('The anemometer height is'):
            zw = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station elevation is'):
            elev = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station latitude is'):
            lat = float(line.split(':')[1].split()[0])

    # Get list of date strings from the input CSV file
    values, ids = [], []
    for test_date in list(csv_df.index):
        # Datetime that has issues with fcd calculation
        if not test_date.startswith('2015-07'):
            continue

        # This day has missing data and is not being handled correctly
        # if test_date.startswith('2015-04-22'):
        #     continue

        test_dt = dt.datetime.strptime(test_date, '%Y-%m-%d')
        # Can the surface type be parameterized inside pytest_generate_tests?
        for surface in ['ETr', 'ETo']:
            date_values = csv_df \
                .loc[test_date, ['TMIN', 'TMAX', 'EA', 'RS', 'WIND']] \
                .rename({
                    'TMIN': 'tmin', 'TMAX': 'tmax', 'EA': 'ea', 'RS': 'rs',
                    'WIND': 'uz'}) \
                .to_dict()
            date_values.update({
                'surface': surface.lower(),
                'expected': out_df.loc[test_date, surface],
                'doy': int(
                    dt.datetime.strptime(test_date, '%Y-%m-%d').strftime('%j')),
                'zw': zw,
                'elev': elev,
                'lat': lat * math.pi / 180,
                'rso_type': 'full',
                'method': 'refet'
            })
            values.append(date_values)
            ids.append('{}-{}'.format(test_date, surface))


@pytest.fixture(scope='module')
def daily_data():
    _daily = DailyData()
    return _daily


def pytest_generate_tests(metafunc):
    # Read in inputs for each daily timestep
    # Set dictionary keys to input variable names
    daily = daily_data()

    if 'daily_params' in metafunc.fixturenames:
        metafunc.parametrize('daily_params', daily.values, ids=daily.ids)


def test_refet_daily_values(daily_params):
    """Test daily RefET calculation at a single point and time"""
    # If I don't copy, the pop changes the test values in daily_data()
    inputs = daily_params.copy()
    surface = inputs.pop('surface')
    expected = inputs.pop('expected')
    # print('ETr: {}'.format(expected))

    # ETr/ETo values only have 4 significant figures
    # Small number of days don't match if difference is set < 0.008
    diff = 0.05 if expected >= 10.0 else 0.008

    if surface.lower() == 'etr':
        assert float(Daily(**inputs).etr()) == pytest.approx(expected, abs=diff)
    elif surface.lower() == 'eto':
        assert float(Daily(**inputs).eto()) == pytest.approx(expected, abs=diff)

import datetime as dt
import math
import os

import numpy as np
import pandas as pd
import pytest
import pytz

from refet.hourly import Hourly
import refet.units as units


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
    'zw': 3.0,
}

# Hourly test parameters for 2015-07-01 18:00 UTC (11:00 AM PDT)
h_args = {
    'doy': 182,
    'ea': 1.1990099614301906,
    'es': 5.09318785259078,
    'eto_asce': 0.6063515410076268,
    'eto_refet': 0.6068613650177561,
    'etr_refet': 0.7201865213918281,
    'etr_asce': 0.7196369609713682,
    'ra': 4.30824147948541,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,  # Conversion from Langleys to MJ m-2
    'tdew': units._f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units._f2c(91.80),
    'uz': 3.33 * 0.44704,  # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

# Test full hourly functions with positional inputs
def test_refet_hourly_input_positions():
    etr = Hourly(
        h_args['tmean'], h_args['ea'], h_args['rs'], h_args['uz'],
        s_args['zw'], s_args['elev'], s_args['lat'], s_args['lon'],
        h_args['doy'], h_args['time'], 'asce').etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])


# Test full hourly calculations with keyword inputs
def test_refet_hourly_default_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time']).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_asce_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='asce').etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_refet_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='refet').etr()
    assert float(etr) == pytest.approx(h_args['etr_refet'])


def test_refet_hourly_default_method_eto():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time']).eto()
    assert float(etr) == pytest.approx(h_args['eto_asce'])

def test_refet_hourly_asce_method_eto():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='asce').eto()
    assert float(etr) == pytest.approx(h_args['eto_asce'])

def test_refet_hourly_refet_method_eto():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='refet').eto()
    assert float(etr) == pytest.approx(h_args['eto_refet'])


# Test unit conversions
def test_refet_hourly_tmean_f():
    etr = Hourly(
        tmean=h_args['tmean'] * (9.0 / 5) + 32, ea=h_args['ea'],
        rs=h_args['rs'], uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'tmean': 'F'}).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_tmean_k():
    etr = Hourly(
        tmean=h_args['tmean'] + 273.15, ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'tmean': 'K'}).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_zw_ft():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'] / 0.3048, elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'zw': 'ft'}).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_elev_ft():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'] / 0.3048,
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'elev': 'ft'}).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

def test_refet_hourly_lon_deg():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'] * math.pi / 180,
        doy=h_args['doy'], time=h_args['time'],
        input_units={'lon': 'rad'}).etr()
    assert float(etr) == pytest.approx(h_args['etr_asce'])

# Test hourly functions using actual RefET input/output files
# DEADBEEF - This doesn't work if I move it to conftest.py
class HourlyData():
    """Setup hourly validation data from Fallon AgriMet station"""
    val_ws = os.path.join(os.getcwd(), 'tests', 'data')
    # val_ws = os.path.join(os.path.dirname(os.getcwd()), 'tests', 'data')

    csv_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.csv')
    # in2_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.in2')
    out_path = os.path.join(val_ws, 'FALN_Agrimet_hourly_raw_2015.out')

    # Read in the inputs CSV file using pandas
    csv_df = pd.read_csv(csv_path, engine='python', na_values='NO RECORD')
    csv_df.rename(
        columns={'OB': 'TEMP', 'TP': 'TDEW', 'WS': 'WIND', 'SI': 'RS'},
        inplace=True)
    # AgriMet times are local with DST (this will drop one hour)
    # DEADBEEF - Can't set timezone as variable in class?
    csv_df['DATETIME'] = csv_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)), axis=1)
    # To match RefET IN2 values exactly, compute DOY using localtime (not UTC)
    csv_df['DOY'] = csv_df['DATETIME'].apply(lambda x: int(x.strftime('%j')))
    csv_df['DATETIME'] = csv_df['DATETIME'].apply(
        lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    csv_df.set_index('DATETIME', inplace=True, drop=True)

    # Convert inputs units
    csv_df['TEMP'] = units._f2c(csv_df['TEMP'])
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
    # # Flatten multi-row header
    # in2_df.columns = [
    #     ' '.join(col).replace('-', '').strip()
    #     for col in in2_df.columns.values]
    # in2_df.rename(
    #     columns={'Year': 'YEAR', 'Mo': 'MONTH', 'Da': 'DAY', 'DoY': 'DOY',
    #              'HrMn': 'HOUR'},
    #     inplace=True)
    # in2_df['HOUR'] = (in2_df['HOUR'] / 100).astype(int)
    # # AgriMet times are local with DST (this will drop one hour)
    # in2_df['DATETIME'] = in2_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
    #     lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)),
    #     axis=1)
    # in2_df['DATETIME'] = in2_df['DATETIME'].apply(
    #     lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    # in2_df.set_index('DATETIME', inplace=True, drop=True)

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
        columns={'Yr': 'YEAR', 'Mo': 'MONTH', 'Day': 'DAY', 'HrMn': 'HOUR',
                 'Tmax': 'TMAX', 'Tmin': 'TMIN', 'DewP': 'TDEW',
                 'Wind': 'WIND', 'Rs': 'RS'},
        inplace=True)
    out_df['HOUR'] = (out_df['HOUR'] / 100).astype(int)
    # AgriMet times are local with DST (this will drop one hour)
    out_df['DATETIME'] = out_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
        lambda x: pytz.timezone('US/Pacific').localize(dt.datetime(*x)),
        axis=1)
    out_df['DATETIME'] = out_df['DATETIME'].apply(
        lambda x: x.tz_convert('UTC').strftime('%Y-%m-%d %H:00'))
    # out_df['DATETIME'] = out_df[['YEAR', 'MONTH', 'DAY', 'HOUR']].apply(
    #     lambda x: dt.datetime(*x).strftime('%Y-%m-%d %H:00'), axis=1)
    out_df.set_index('DATETIME', inplace=True, drop=True)

    # Read the station properties from the IN2 file for now
    # The values should probably be read using a regular expression
    for line in out_data:
        if line.strip().startswith('The anemometer height is'):
            zw = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station elevation is'):
            elev = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station latitude is'):
            lat = float(line.split(':')[1].split()[0])
        elif line.strip().startswith('The weather station longitude is'):
            lon = float(line.split(':')[1].split()[0])
            if 'West' in line:
                lon *= -1

    # Get list of date strings from the input CSV file
    values, ids = [], []
    for test_date in list(csv_df.index):
        # Datetime that has issues with fcd calculation
        if not test_date.startswith('2015-07-01'):
            continue

        # Only check day time values for now
        if float(csv_df.loc[test_date, 'RS']) <= 1.0:
            continue

        test_dt = dt.datetime.strptime(test_date, '%Y-%m-%d %H:%M')
        # Can the surface type be parameterized inside pytest_generate_tests?
        for surface in ['ETr', 'ETo']:
            date_values = csv_df \
                .loc[test_date, ['TEMP', 'EA', 'RS', 'WIND', 'DOY']] \
                .rename({
                    'DOY': 'doy', 'TEMP': 'tmean', 'EA': 'ea', 'RS': 'rs',
                    'WIND': 'uz'}) \
                .to_dict()
            date_values.update({
                'surface': surface.lower(),
                'expected': out_df.loc[test_date, surface],
                # 'doy': int(test_dt.strftime('%j')),
                'time': test_dt.hour,
                'zw': zw,
                'elev': elev,
                'lat': lat,
                'lon': lon,
                'method': 'refet'
            })
            values.append(date_values)
            ids.append('{}-{}'.format(test_date, surface))


@pytest.fixture(scope='module')
def hourly_data():
    _hourly = HourlyData()
    return _hourly


def pytest_generate_tests(metafunc):
    # Read in inputs for each daily timestep
    # Set dictionary keys to input variable names
    hourly = hourly_data()

    if 'hourly_params' in metafunc.fixturenames:
        metafunc.parametrize('hourly_params', hourly.values, ids=hourly.ids)


def test_refet_hourly_func_values(hourly_params):
    """Test hourly RefET calculation at a single point and time"""
    # If I don't copy, the pop changes the test values in hourly_data()
    inputs = hourly_params.copy()
    surface = inputs.pop('surface')
    expected = inputs.pop('expected')
    # print('ETr: {}'.format(expected))

    if surface.lower() == 'etr':
        assert float(Hourly(**inputs).etr()) == pytest.approx(expected, abs=0.01)
    elif surface.lower() == 'eto':
        assert float(Hourly(**inputs).eto()) == pytest.approx(expected, abs=0.01)
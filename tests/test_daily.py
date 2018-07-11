import math

import pytest

from refet.daily import Daily
import refet.units as units


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': 39.4575,
    'lon': -118.77388,
    'zw': 3.0,
}

# Daily test parameters for 2015-07-01
# Ea was computed from q for both asce and refet methods
d_args = {
    'doy': 182,
    'ea': 1.2206674169951346,
    'eto_asce': 7.9422320475712835,
    'eto_refet': 7.9422320475712835,
    'etr_asce': 10.626087665395694,
    'etr_refet': 10.571314344056955,
    'etr_rso_simple': 10.628137858930051,
    'rs': 674.07 * 0.041868,  # Convert Rs from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'uz': 4.80 * 0.44704,  # Conversion wind speed from mph to m s-1
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


# Test rso parameters
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


@pytest.mark.parametrize(
    'method, expected',
    [['asce', d_args['etr_asce']],
     ['refet', d_args['etr_refet']]])
def test_refet_daily_method(method, expected):
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        method=method).etr()
    assert float(etr) == pytest.approx(expected)


@pytest.mark.parametrize(
    'surface, expected',
    [['etr', d_args['etr_refet']],
     ['alfalfa', d_args['etr_refet']],
     ['tall', d_args['etr_refet']],
     ['eto', d_args['eto_refet']],
     ['grass', d_args['eto_refet']],
     ['short', d_args['eto_refet']]])
def test_refet_daily_etsz(surface, expected):
    etsz = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'], method='refet').etsz(surface)
    assert float(etsz) == pytest.approx(expected)


# Need to test other input structures (i.e. 1d and 2d arrays)
# def test_refet_daily_1d():
#     assert True

# def test_refet_daily_2d():
#     assert True


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


def test_refet_daily_ea_pa():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'] * 1000,
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'ea': 'Pa'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_rs_langleys():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'] / 0.041868, uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'rs': 'Langleys'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_rs_wm2():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'] / 0.0864, uz=d_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], doy=d_args['doy'],
        input_units={'rs': 'W m-2'}).etr()
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


def test_refet_daily_lat_deg():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], doy=d_args['doy'],
        input_units={'lat': 'deg'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])


def test_refet_daily_lat_rad():
    etr = Daily(
        tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
        rs=d_args['rs'], uz=d_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'] * math.pi / 180, doy=d_args['doy'],
        input_units={'lat': 'rad'}).etr()
    assert float(etr) == pytest.approx(d_args['etr_asce'])

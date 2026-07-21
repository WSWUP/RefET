import math

import pytest

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
    'q': 0.008536365803069757,
    'ra': 4.30824147948541,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,  # Convert Rs from Langleys to MJ m-2
    'tdew': units.f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units.f2c(91.80),
    'uz': 3.33 * 0.44704,  # Conversion wind speed from mph to m s-1
    'u2': 1.3709275319197722,
}

# Test full hourly functions with positional inputs
def test_refet_hourly_input_positions():
    etr = Hourly(
        h_args['tmean'], h_args['rs'], h_args['uz'],
        s_args['zw'], s_args['elev'], s_args['lat'], s_args['lon'],
        h_args['doy'], h_args['time'], h_args['ea'], None, 'asce',
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


# Test full hourly calculations with keyword inputs
def test_refet_hourly_default_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'],
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_asce_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='asce',
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_refet_method_etr():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='refet',
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_refet'])


def test_refet_hourly_default_method_eto():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'],
    ).eto()
    assert float(etr[0]) == pytest.approx(h_args['eto_asce'])


@pytest.mark.parametrize(
    'method, expected',
    [
        ['asce', h_args['eto_asce']],
        ['refet', h_args['eto_refet']],
    ]
)
def test_refet_hourly_asce_method_eto(method, expected):
    eto = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method=method,
    ).eto()
    assert float(eto[0]) == pytest.approx(expected)


@pytest.mark.parametrize(
    'surface, expected',
    [
        ['etr', h_args['etr_refet']],
        ['alfalfa', h_args['etr_refet']],
        ['tall', h_args['etr_refet']],
        ['eto', h_args['eto_refet']],
        ['grass', h_args['eto_refet']],
        ['short', h_args['eto_refet']],
    ]
)
def test_refet_daily_etsz(surface, expected):
    etsz = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], method='refet',
    ).etsz(surface)
    assert float(etsz[0]) == pytest.approx(expected)


# Test unit conversions
def test_refet_hourly_tmean_f():
    etr = Hourly(
        tmean=h_args['tmean'] * (9.0 / 5) + 32, ea=h_args['ea'],
        rs=h_args['rs'], uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'tmean': 'F'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_tmean_k():
    etr = Hourly(
        tmean=h_args['tmean'] + 273.15, ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'tmean': 'K'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_ea_pa():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'] * 1000, rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'ea': 'Pa'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_rs_langleys():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'] / 0.041868,
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'rs': 'Langleys'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_rs_wm2():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'] / 0.0036,
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'rs': 'W m-2'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_zw_ft():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'] / 0.3048, elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'zw': 'ft'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_elev_ft():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'] / 0.3048,
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'elev': 'ft'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_lon_deg():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'],
        doy=h_args['doy'], time=h_args['time'],
        input_units={'lon': 'deg'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_lon_rad():
    etr = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'] * math.pi / 180,
        doy=h_args['doy'], time=h_args['time'],
        input_units={'lon': 'rad'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_tdew():
    """Check that Ea can be computed from tdew"""
    etr = Hourly(
        tmean=h_args['tmean'], tdew=h_args['tdew'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'],
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


def test_refet_hourly_ea_tdew_exception():
    """Check that a ValueError is raised if ea and tdew are not set"""
    with pytest.raises(ValueError):
        Hourly(
            tmean=h_args['tmean'], rs=h_args['rs'],
            uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
            lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
            time=h_args['time'],
        )


def test_refet_hourly_tdew_f():
    """Check that Tdew is converted to C before computing Ea"""
    etr = Hourly(
        tmean=h_args['tmean'], rs=h_args['rs'], tdew=units.c2f(h_args['tdew']),
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'], input_units={'tdew': 'F'},
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])


# @pytest.mark.parametrize(
#     'tmean, ea, rs, uz, zw, elev, lat, lon, doy, time, method, expected',
#     [
#         [30.56, 315.0, 704.5 * 0.0036, 3.5, 2.0, 8.0, 39.0, -120.0, 183, 20.0, 'asce', 1],
#         # [30.56, 315.0, 704.5 * 0.0036, 3.5, 2.0, 8.0, 15.2, 145.7, 183, 0, 'asce', 1],
#     ]
# )
# def test_refet_hourly_values(tmean, ea, rs, uz, zw, elev, lat, lon, doy, time, method, expected):
#     eto = Hourly(
#         tmean=tmean, ea=ea, rs=rs, uz=uz, zw=zw, elev=elev, lat=lat, lon=lon,
#         doy=doy, time=time, method=method,
#         input_units={'lon': 'deg', 'lat': 'deg', 'rs': 'w m-2'},
#     ).eto()
#     assert float(eto[0]) == pytest.approx(expected)


def test_refet_hourly_doy_is_array():
    """doy must survive as an ndarray (it was previously overwritten)"""
    import numpy as np
    h = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=[h_args['doy']],
        time=h_args['time'],
    )
    assert isinstance(h.doy, np.ndarray)


def test_refet_hourly_vpd_clamped():
    """VPD is clamped at zero when Tdew exceeds Tmean, matching Daily"""
    h = Hourly(
        tmean=10.0, tdew=12.0, rs=0.0, uz=h_args['uz'], zw=s_args['zw'],
        elev=s_args['elev'], lat=s_args['lat'], lon=s_args['lon'],
        doy=h_args['doy'], time=8.0,
    )
    assert float(h.vpd[0]) == 0.0



def test_refet_hourly_results():
    """results() returns the constructor's intermediate terms by name"""
    h = Hourly(
        tmean=h_args['tmean'], ea=h_args['ea'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'])
    r = h.results()
    assert set(r.keys()) == {
        'tmean', 'ea', 'es', 'es_slope', 'vpd',
        'rs', 'pair', 'psy', 'ra', 'rso', 'fcd', 'rnl', 'rn', 'u2'}
    assert float(r['u2'][0]) == pytest.approx(h_args['u2'])
    assert r['ra'] is h.ra


def test_refet_hourly_q():
    """Ea can be computed from specific humidity"""
    etr = Hourly(
        tmean=h_args['tmean'], q=h_args['q'], rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'],
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'], abs=0.001)


def test_refet_hourly_rh():
    """Ea can be computed from relative humidity (ASCE Table 4)"""
    rh = 100 * h_args['ea'] / h_args['es']
    etr = Hourly(
        tmean=h_args['tmean'], rh=rh, rs=h_args['rs'],
        uz=h_args['uz'], zw=s_args['zw'], elev=s_args['elev'],
        lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
        time=h_args['time'],
    ).etr()
    assert float(etr[0]) == pytest.approx(h_args['etr_asce'])

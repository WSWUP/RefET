import math

import numpy as np
import pytest

import src.refet.calcs as calcs
import src.refet.units as units


# Eventually move to conftest.py or a separate file
# Fallon AgriMet site parameters
s_args = {
    'elev': 1208.5,
    'lat': units._deg2rad(39.4575),
    'lon': units._deg2rad(-118.77388),
    'pair': 87.81876435813037,
    'pair_asce': 87.80710537212929,
    'zw': 3.0,
}
# Daily test parameters for 2015-07-01
d_args = {
    'doy': 182,
    'delta': 0.40352881013673136,
    'delta_asce': 0.4029517192078854,
    'doy_frac': 3.132985550429273,
    'dr': 0.9670012223491632,
    'ea': 1.2206674169951346,           # Computed from Tdew
    'es': 4.6747236227258835,
    'es_slope': 0.23489129849801055,
    'es_slope_asce': 0.23488581814172638,
    'eto': 7.942481120179387,
    'etr': 10.571560006380153,
    'fcd': 0.8569860867772078,
    'omega_s': 1.9298904620748385,
    'q': 0.008691370735727117,          # Computed from Ea from Tdew
    'q_asce': 0.008692530868140688,     # Computed from Ea from Tdew
    'ra': 41.67610845067083,
    'ra_asce': 41.64824567735701,
    'rn': 15.174377350374275,
    'rnl': 6.556533974825727,
    'rs': 674.07 * 0.041868,            # Conversion from Langleys to MJ m-2
    'rso': 31.565939444861765,
    'rso_simple': 32.26439287925584,
    'sc': -0.05874166519510547,
    'tdew': units._f2c(49.84),
    'tmin': units._f2c(66.65),
    'tmax': units._f2c(102.80),
    # 'tmean': f2c(84.725),
    'w': 17.107650595384076,
    'uz': 4.80 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.976111757722194,
}
# Hourly test parameters for 2015-07-01 18:00 UTC (11:00 AM PDT)
h_args = {
    'doy': 182,
    'ea': 1.1990099614301906,           # Computed from Tdew
    'es': 5.09318785259078,
    'es_slope': 0.2852040944519568,
    'es_slope_asce': 0.28519744022482574,
    'eto': 0.6065255163817055,
    'etr': 0.7201865213918281,
    'fcd': 0.6816142001345745,
    'fcd_asce': 0.6816142001345745,
    # 'omega': -0.3866777826605525,     # Computed from time_mid
    'omega': -0.5175774765601271,       # Computed from time to match IN2
    'q': 0.008536365803069757,          # Computed from Ea from Tdew
    'q_asce': 0.00853750513849305,      # Computed from Ea from Tdew
    'ra': 4.30824147948541,
    'ra_asce': 4.30635461533285,
    'rn': 1.7427193535884922,
    'rnl': 0.22897874401150786,
    'rs': 61.16 * 0.041868,             # Conversion from Langleys to MJ m-2
    'rso': 3.350936122776373,
    'rso_asce': 3.350851392549374,
    'rso_simple': 3.33531130617322,
    # 'solar_time': -1.4770003318617722,  # Computed from time_mid
    'solar_time': -1.9770003318617722,  # Computed from time to match IN2
    'tdew': units._f2c(49.36),
    'time': 18.0,
    'time_mid': 18.5,
    'tmean': units._f2c(91.80),
    'uz': 3.33 * 0.44704,               # Conversion from mph to m s-1
    'u2': 1.3709275319197722,
}

# # Playing around with passing a dictionary of function keyword arguments
# daily_args = {
#     k: v for k, v in daily.items()
#     if k in ['tmin', 'tmax', 'ea', 'rs', 'uz', 'zw', 'doy']}
# daily_args.update({'ref_type':'etr', 'rso_type': 'full', 'rso': None})
# hourly_args = {
#     k: v for k, v in daily.items()
#     if k in ['tmean', 'ea', 'rs', 'uz', 'zw', 'doy', 'time']}
# hourly_args.update({'ref_type':'etr'})


## Test ancillary functions with positional inputs
def test_air_pressure_default(elev=s_args['elev'], pair=s_args['pair_asce']):
    assert float(calcs._air_pressure(elev)) == pytest.approx(pair)

def test_air_pressure_asce(elev=s_args['elev'], pair=s_args['pair_asce']):
    assert float(calcs._air_pressure(elev, method='asce')) == pytest.approx(pair)

def test_air_pressure_refet(elev=s_args['elev'], pair=s_args['pair']):
    assert float(calcs._air_pressure(elev, method='refet')) == pytest.approx(pair)


@pytest.mark.parametrize(
    'tdew, ea',
    [
        [d_args['tdew'], d_args['ea']],
        [h_args['tmean'], h_args['es']],
        [h_args['tdew'], h_args['ea']],
    ]
)
def test_sat_vapor_pressure(tdew, ea):
    assert float(calcs._sat_vapor_pressure(tdew)) == pytest.approx(ea)


@pytest.mark.parametrize(
    'ea, pair, q',
    [
        [d_args['ea'], s_args['pair'], d_args['q']],
        [d_args['ea'], s_args['pair_asce'], d_args['q_asce']],
        [h_args['ea'], s_args['pair'], h_args['q']],
        [h_args['ea'], s_args['pair_asce'], h_args['q_asce']],
    ]
)
def test_specific_humidity_daily(ea, pair, q):
    assert float(calcs._specific_humidity(ea, pair)) == pytest.approx(q)


def test_actual_vapor_pressure(q=d_args['q'], pair=s_args['pair'],
                               ea=d_args['ea']):
    assert float(calcs._actual_vapor_pressure(q, pair)) == pytest.approx(ea)


def test_vpd(es=d_args['es'], ea=d_args['ea']):
    assert float(calcs._vpd(es, ea)) == pytest.approx(float(es-ea))
    # Check that negative VPD's are set to 0
    assert float(calcs._vpd(es, es+1)) == pytest.approx(0)


def test_es_slope_default(tmin=d_args['tmin'], tmax=d_args['tmax'],
                          es_slope=d_args['es_slope_asce']):
    output = calcs._es_slope(0.5 * (tmin + tmax))
    assert float(output) == pytest.approx(es_slope)


def test_es_slope_asce(tmin=d_args['tmin'], tmax=d_args['tmax'],
                       es_slope=d_args['es_slope_asce']):
    output = calcs._es_slope(0.5 * (tmin + tmax), method='asce')
    assert float(output) == pytest.approx(es_slope)


def test_es_slope_refet(tmin=d_args['tmin'], tmax=d_args['tmax'],
                        es_slope=d_args['es_slope']):
    output = calcs._es_slope(0.5 * (tmin + tmax), method='refet')
    assert float(output) == pytest.approx(es_slope)


def test_es_slope_hourly_asce(tmean=h_args['tmean'], es_slope=h_args['es_slope_asce']):
    assert float(calcs._es_slope(tmean, method='asce')) == pytest.approx(es_slope)


def test_es_slope_hourly_refet(tmean=h_args['tmean'], es_slope=h_args['es_slope']):
    assert float(calcs._es_slope(tmean, method='refet')) == pytest.approx(es_slope)


def test_precipitable_water(pair=s_args['pair'], ea=d_args['ea'], w=d_args['w']):
    assert float(calcs._precipitable_water(pair, ea)) == pytest.approx(w)


def test_doy_fraction(doy=d_args['doy'], expected=d_args['doy_frac']):
    assert float(calcs._doy_fraction(doy)) == pytest.approx(expected)


def test_delta_default(doy=d_args['doy'], delta=d_args['delta_asce']):
    assert float(calcs._delta(doy)) == pytest.approx(delta)


def test_delta_asce(doy=d_args['doy'], delta=d_args['delta_asce']):
    assert float(calcs._delta(doy, method='asce')) == pytest.approx(delta)


def test_delta_refet(doy=d_args['doy'], delta=d_args['delta']):
    assert float(calcs._delta(doy, method='refet')) == pytest.approx(delta)


def test_dr(doy=d_args['doy'], dr=d_args['dr']):
    assert float(calcs._dr(doy)) == pytest.approx(dr)


def test_seasonal_correction(doy=d_args['doy'], sc=d_args['sc']):
    assert float(calcs._seasonal_correction(doy)) == pytest.approx(sc)


def test_solar_time_rad(lon=s_args['lon'], time_mid=h_args['time'],
                        sc=d_args['sc'], expected=h_args['solar_time']):
    assert float(calcs._solar_time_rad(lon, time_mid, sc)) == pytest.approx(expected)


def test_omega(solar_time=h_args['solar_time'], omega=h_args['omega']):
    assert float(calcs._omega(solar_time)) == pytest.approx(omega)


@pytest.mark.parametrize(
    'x, x_min, x_max, expected',
    [
        [1.2, 1.2, 1.5, 1.2],
        [1.1, 1.2, 1.5, 1.4],
        [1.6, 1.2, 1.5, 1.3],
        [2.0, 1.2, 1.5, 1.4],
    ]
)
def test_wrap(x, x_min, x_max, expected):
    value = calcs._wrap(x, x_min, x_max)
    assert pytest.approx(float(value)) == expected


def test_omega_sunset(lat=s_args['lat'], delta=d_args['delta'],
                      omega_s=d_args['omega_s']):
    assert float(calcs._omega_sunset(lat, delta)) == pytest.approx(omega_s)


def test_omega_sunset_high_lat():
    # Test that omega sunset is limited to [0, pi]
    # This occurs for angles > ~65 degrees in the summer and ~75 in the winter
    lat = 70 * (math.pi / 180)
    assert float(calcs._omega_sunset(lat, calcs._delta(182))) == pytest.approx(math.pi)
    assert float(calcs._omega_sunset(-lat, calcs._delta(182))) == pytest.approx(0)
    assert float(calcs._omega_sunset(lat, calcs._delta(1))) == pytest.approx(0)


def test_ra_daily_default(lat=s_args['lat'], doy=d_args['doy'], ra=d_args['ra_asce']):
    assert float(calcs._ra_daily(lat, doy)) == pytest.approx(ra)


def test_ra_daily_asce(lat=s_args['lat'], doy=d_args['doy'], ra=d_args['ra_asce']):
    assert float(calcs._ra_daily(lat, doy, method='asce')) == pytest.approx(ra)


def test_ra_daily_refet(lat=s_args['lat'], doy=d_args['doy'], ra=d_args['ra']):
    assert float(calcs._ra_daily(lat, doy, method='refet')) == pytest.approx(ra)


def test_ra_daily_zero():
    # Ra can go to zero for winter DOY and/or high latitudes
    assert float(calcs._ra_daily(68 * (math.pi / 180), 1)) == 0
    assert float(calcs._ra_daily(80 * (math.pi / 180), 55)) == 0


def test_ra_hourly_default(lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
                           time=h_args['time_mid'], ra=h_args['ra_asce']):
    assert float(calcs._ra_hourly(lat, lon, doy, time)) == pytest.approx(ra)

def test_ra_hourly_asce(lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
                        time=h_args['time_mid'], ra=h_args['ra_asce']):
    output = calcs._ra_hourly(lat, lon, doy, time, method='asce')
    assert float(output) == pytest.approx(ra)

def test_ra_hourly_refet(lat=s_args['lat'], lon=s_args['lon'], doy=h_args['doy'],
                         time=h_args['time_mid'], ra=h_args['ra']):
    output = calcs._ra_hourly(lat, lon, doy, time, method='refet')
    assert float(output) == pytest.approx(ra)


def test_rso_daily(ra=d_args['ra'], ea=d_args['ea'], pair=s_args['pair'],
                   doy=d_args['doy'], lat=s_args['lat'], rso=d_args['rso']):
    assert float(calcs._rso_daily(ra, ea, pair, doy, lat)) == pytest.approx(rso)


def test_rso_daily_ra_zero(ea=d_args['ea'], pair=s_args['pair']):
    # Rso can go to zero for winter DOY and/or high latitudes when Ra is zero
    output = calcs._rso_daily(
        calcs._ra_daily(80 * math.pi / 180, 1), ea, pair, 1, 80 * math.pi / 180
    )
    assert float(output) == 0


def test_rso_hourly_default(ra=h_args['ra'], ea=h_args['ea'], pair=s_args['pair'],
                            doy=h_args['doy'], time=h_args['time_mid'],
                            lat=s_args['lat'], lon=s_args['lon'],
                            rso=h_args['rso_asce']):
    output = calcs._rso_hourly(ra, ea, pair, doy, time, lat, lon)
    assert float(output) == pytest.approx(rso)


def test_rso_hourly_asce(ra=h_args['ra'], ea=h_args['ea'], pair=s_args['pair'],
                         doy=h_args['doy'], time=h_args['time_mid'],
                         lat=s_args['lat'], lon=s_args['lon'],
                         rso=h_args['rso_asce']):
    output = calcs._rso_hourly(ra, ea, pair, doy, time, lat, lon,method='asce')
    assert float(output) == pytest.approx(rso)


def test_rso_hourly_refet(ra=h_args['ra'], ea=h_args['ea'], pair=s_args['pair'],
                          doy=h_args['doy'], time=h_args['time_mid'],
                          lat=s_args['lat'], lon=s_args['lon'],
                          rso=h_args['rso']):
    output = calcs._rso_hourly(ra, ea, pair, doy, time, lat, lon, method='refet')
    assert float(output) == pytest.approx(rso)


@pytest.mark.parametrize(
    'ra, elev, rso',
    [
        [d_args['ra'], s_args['elev'], d_args['rso_simple']],
        [h_args['ra'], s_args['elev'], h_args['rso_simple']],
    ]
)
def test_rso_simple(ra, elev, rso):
    assert float(calcs._rso_simple(ra, elev)) == pytest.approx(rso)


def test_fcd_daily(rs=d_args['rs'], rso=d_args['rso'], fcd=d_args['fcd']):
    assert float(calcs._fcd_daily(rs, rso)) == pytest.approx(fcd)


def test_fcd_daily_rso_zero():
    # Check that fcd returns 1 when Rso is zero
    assert float(calcs._fcd_daily(1, 0)) == pytest.approx(1.0)
    assert float(calcs._fcd_daily(0, 0)) == pytest.approx(1.0)


def test_fcd_hourly_default(rs=h_args['rs'], rso=h_args['rso'],
                            doy=h_args['doy'], time=h_args['time_mid'],
                            lat=s_args['lat'], lon=s_args['lon'],
                            fcd=h_args['fcd_asce']):
    output = calcs._fcd_hourly(rs, rso, doy, time, lat, lon)
    assert float(output) == pytest.approx(fcd)


def test_fcd_hourly_asce(rs=h_args['rs'], rso=h_args['rso'],
                         doy=h_args['doy'], time=h_args['time_mid'],
                         lat=s_args['lat'], lon=s_args['lon'],
                         fcd=h_args['fcd_asce']):
    output = calcs._fcd_hourly(rs, rso, doy, time, lat, lon, method='asce')
    assert float(output) == pytest.approx(fcd)


def test_fcd_hourly_refet(rs=h_args['rs'], rso=h_args['rso'],
                          doy=h_args['doy'], time=h_args['time_mid'],
                          lat=s_args['lat'], lon=s_args['lon'],
                          fcd=h_args['fcd']):
    output = calcs._fcd_hourly(rs, rso, doy, time, lat, lon, method='refet')
    assert float(output) == pytest.approx(fcd)


def test_fcd_hourly_night(rs=h_args['rs'], rso=h_args['rso'],
                          doy=h_args['doy'], time=6,
                          lat=s_args['lat'], lon=s_args['lon'], fcd=1):
    # For now, check that nighttime fcd values are set to 1
    output = calcs._fcd_hourly(rs, rso, doy, time, lat, lon, method='refet')
    assert float(output) == pytest.approx(fcd)


# Test function when rs/rso have different shapes than lat/lon/doy/time
# def test_fcd_hourly_broadcasting():
#     assert False


def test_rnl_daily(tmin=d_args['tmin'], tmax=d_args['tmax'], ea=d_args['ea'],
                   fcd=d_args['fcd'], rnl=d_args['rnl']):
    assert float(calcs._rnl_daily(tmax, tmin, ea, fcd)) == pytest.approx(rnl)


def test_rnl_hourly(tmean=h_args['tmean'], ea=h_args['ea'],
                    fcd=h_args['fcd'], rnl=h_args['rnl']):
    assert float(calcs._rnl_hourly(tmean, ea, fcd)) == pytest.approx(rnl)


def test_rn_daily(rs=d_args['rs'], rnl=d_args['rnl'], rn=d_args['rn']):
    assert float(calcs._rn_daily(rs, rnl)) == pytest.approx(rn)


def test_rn_hourly(rs=h_args['rs'], rnl=h_args['rnl'], rn=h_args['rn']):
    assert float(calcs._rn_hourly(rs, rnl)) == pytest.approx(rn)


@pytest.mark.parametrize(
    'uz, zw, u2',
    [
        [d_args['uz'], s_args['zw'], d_args['u2']],
        [h_args['uz'], s_args['zw'], h_args['u2']],
    ]
)
def test_wind_height_adjust(uz, zw, u2):
    assert float(calcs._wind_height_adjust(uz, zw)) == pytest.approx(u2)


def test_wind_height_adjust_2m(uz=2.5, zw=2.0, u2=2.5):
    assert float(calcs._wind_height_adjust(uz, zw)) == pytest.approx(u2, abs=0.001)


@pytest.mark.parametrize(
    'rn, g, tmean, u2, es, ea, es_slope, pair, cn, cd, etsz',
    [
        [d_args['rn'], 0, 0.5 * (d_args['tmin'] + d_args['tmax']), d_args['u2'],
         d_args['es'], d_args['ea'], d_args['es_slope'], s_args['pair'],
         900, 0.34, d_args['eto']],
        [d_args['rn'], 0, 0.5 * (d_args['tmin'] + d_args['tmax']), d_args['u2'],
         d_args['es'], d_args['ea'], d_args['es_slope'], s_args['pair'],
         1600, 0.38, d_args['etr']],
        [h_args['rn'], 0.1 * h_args['rn'], h_args['tmean'], h_args['u2'],
         h_args['es'], h_args['ea'], h_args['es_slope'], s_args['pair'],
         37.0, 0.24, h_args['eto']],
        [h_args['rn'], 0.04 * h_args['rn'], h_args['tmean'], h_args['u2'],
         h_args['es'], h_args['ea'], h_args['es_slope'], s_args['pair'],
         66.0, 0.25, h_args['etr']],
    ]
)
def test_etsz(rn, g, tmean, u2, es, ea, es_slope, pair, cn, cd, etsz):
    output = calcs._etsz(rn, g, tmean, u2, es - ea, es_slope, 0.000665 * pair, cn, cd)
    assert float(output) == pytest.approx(etsz, abs=0.001)

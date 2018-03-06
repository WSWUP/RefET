import math

import numpy as np


def _air_pressure(elev, method='refet'):
    """Mean atmospheric pressure at station elevation (Eqs. 3 & 34)

    Parameters
    ----------
    elev : scalar or array_like of shape(M, )
        Elevation [m].
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.

    Returns
    -------
    pair : ndarray
        Air pressure [kPa].

    Notes
    -----
    The current calculation in Ref-ET:
        101.3 * (((293 - 0.0065 * elev) / 293) ** (9.8 / (0.0065 * 286.9)))
    Equation 3 in ASCE-EWRI 2005:
        101.3 * (((293 - 0.0065 * elev) / 293) ** 5.26)
    Per Dr. Allen, the calculation with full precision:
        101.3 * (((293.15 - 0.0065 * elev) / 293.15) ** (9.80665 / (0.0065 * 286.9)))

    """
    pair = np.array(elev, copy=True, ndmin=1).astype(np.float64)
    pair *= -0.0065
    if method == 'asce':
        pair += 293
        pair /= 293
        np.power(pair, 5.26, out=pair)
    else:
        pair += 293
        pair /= 293
        np.power(pair, 9.8 / (0.0065 * 286.9), out=pair)
    # np.power(pair, 5.26, out=pair)
    pair *= 101.3
    return pair


def _sat_vapor_pressure(temperature):
    """Saturation vapor pressure from temperature (Eq. 7)

    Parameters
    ----------
    temperature : scalar or array_like of shape(M, )
        Air temperature [C].

    Returns
    -------
    e : ndarray
        Saturation vapor pressure [kPa].

    Notes
    -----
    0.6108 * exp(17.27 * temperature / (temperature + 237.3))

    """
    e = np.array(temperature, copy=True, ndmin=1).astype(np.float64)
    e += 237.3
    np.reciprocal(e, out=e)
    e *= temperature
    e *= 17.27
    np.exp(e, out=e)
    e *= 0.6108
    return e



def _es_slope(tmean, method='refet'):
    """Slope of the saturation vapor pressure-temperature curve (Eq. 5)

    Parameters
    ----------
    tmean : ndarray
        Mean air temperature [C].
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.

    Returns
    -------
    ndarray

    Notes
    -----
    4098 * 0.6108 * exp(17.27 * T / (T + 237.3)) / ((T + 237.3) ** 2))

    """
    if method == 'refet':
        es_slope = (
            4098.0 * _sat_vapor_pressure(tmean) / np.power(tmean + 237.3, 2))
    elif method == 'asce':
        es_slope = (
            2503.0 * np.exp(17.27 * tmean / (tmean + 237.3)) /
            np.power(tmean + 237.3, 2))
    return es_slope


def _actual_vapor_pressure(q, pair):
    """"Actual vapor pressure from specific humidity

    Parameters
    ----------
    q : scalar or array_like of shape(M, )
        Specific humidity [kg/kg].
    pair : scalar or array_like of shape(M, )
        Air pressure [kPa].

    Returns
    -------
    ea : ndarray
        Actual vapor pressure [kPa].

    Notes
    -----
    q * pair / (0.622 + 0.378 * q)

    """
    ea = np.array(q, copy=True, ndmin=1).astype(np.float64)
    ea *= 0.378
    ea += 0.622
    np.reciprocal(ea, out=ea)
    ea *= pair
    ea *= q
    return ea


def _specific_humidity(ea, pair):
    """"Specific humidity from actual vapor pressure

    Parameters
    ----------
    ea : scalar or array_like of shape(M, )
        Specific humidity [kPa].
    pair : scalar or array_like of shape(M, )
        Air pressure [kPa].

    Returns
    -------
    q : ndarray
        Specific humidity [kg/kg].

    Notes
    -----
    0.622 * ea / (pair - 0.378 * ea)

    """
    q = np.array(ea, copy=True, ndmin=1).astype(np.float64)
    q *= -0.378
    q += pair
    np.reciprocal(q, out=q)
    q *= ea
    q *= 0.622
    return q


def _vpd(es, ea):
    """Vapor pressure deficit

    Parameters
    ----------
    es : scalar or array_like of shape(M, )
        Saturated vapor pressure [kPa].
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].

    Returns
    -------
    ndarray
        Vapor pressure deficit [kPa].

    """

    return np.maximum(es - ea, 0)


def _precipitable_water(pair, ea):
    """Precipitable water in the atmosphere (Eq. D.3)

    Parameters
    ----------
    pair : scalar or array_like of shape(M, )
        Air pressure [kPa].
    ea : scalar or array_like of shape(M, )
        Vapor pressure [kPa].

    Returns
    -------
    ndarray
        Precipitable water [mm].

    """
    return pair * 0.14 * ea + 2.1


def _doy_fraction(doy):
    """Fraction of the DOY in the year (Eq. 50)

    Parameters
    ----------
    doy : scalar or array_like of shape(M, )
        Day of year.

    Returns
    -------
    array_like
        DOY fraction [radians].

    """
    return doy * (2.0 * math.pi / 365)


def _delta(doy, method='refet'):
    """Earth declination (Eq. 51)

    Parameters
    ----------
    doy : scalar or array_like of shape(M, )
        Day of year.
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.

    Returns
    -------
    ndarray
        Earth declination [radians].

    Notes
    -----
    Original equation in Duffie & Beckman (1980) (with conversions to radians):
        23.45 * (pi / 180) * sin(2 * pi * (doy + 284) / 365)
    Equation 24 in ASCE-EWRI (2005):
        0.409 * sin((2 * pi * doy / 365) - 1.39)

    """
    if method == 'asce':
        return 0.409 * np.sin(_doy_fraction(doy) - 1.39)
    else:
        return 23.45 * (math.pi / 180) * np.sin(2 * math.pi * (doy + 284) / 365)


def _dr(doy):
    """Inverse square of the Earth-Sun Distance (Eq. 50)

    Parameters
    ----------
    doy : scalar or array_like of shape(M, )
        Day of year.

    Returns
    -------
    ndarray

    Notes
    -----
    This function returns 1 / d^2, not d, for direct use in radiance to
      TOA reflectance calculation
    pi * L * d^2 / (ESUN * cos(theta)) -> pi * L / (ESUN * cos(theta) * d)

    """
    return 1.0 + 0.033 * np.cos(_doy_fraction(doy))


def _seasonal_correction(doy):
    """Seasonal correction for solar time (Eqs. 57 & 58)

    Parameters
    ----------
    doy : scalar or array_like of shape(M, )
        Day of year.

    Returns
    ------
    ndarray
        Seasonal correction [hour]

    """
    b = 2 * math.pi * (doy - 81.) / 364.
    return 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.0250 * np.sin(b)


def _solar_time_rad(lon, time_mid, sc):
    """Solar time (i.e. noon is 0) (Eq. 55)

    Parameters
    ----------
    lon : scalar or array_like of shape(M, )
        Longitude [radians].
    time_mid : scalar or array_like of shape(M, )
        UTC time at midpoint of period [hours].
    sc : scalar or array_like of shape(M, )
        Seasonal correction [hours].

    Returns
    -------
    ndarray
        Solar time [hours].

    Notes
    -----
    This function could be integrated into the _omega() function since they are
    always called together (i.e. _omega(_solar_time_rad()).  It was built
    independently from _omega to eventually support having a separate
    solar_time functions for longitude in degrees.

    """
    return time_mid + (lon * 24 / (2 * math.pi)) + sc - 12


def _omega(solar_time):
    """Solar hour angle (Eq. 55)

    Parameters
    ----------
    solar_time : scalar or array_like of shape(M, )
        Solar time (i.e. noon is 0) [hours].

    Returns
    -------
    omega : ndarray
        Hour angle [radians].

    """
    omega = (2 * math.pi / 24.0) * solar_time

    # Need to adjust omega so that the values go from -pi to pi
    # Values outside this range are wrapped (i.e. -3*pi/2 -> pi/2)
    omega = _wrap(omega, -math.pi, math.pi)
    return omega


def _wrap(x, x_min, x_max):
    """Wrap floating point values into range

    Parameters
    ----------
    x : ndarray
        Values to wrap.
    x_min : float
        Minimum value in output range.
    x_max : float
        Maximum value in output range.

    Returns
    -------
    ndarray

    """
    return np.mod((x - x_min), (x_max - x_min)) + x_min


def _omega_sunset(lat, delta):
    """Sunset hour angle (Eq. 59)

    Parameters
    ----------
    lat : scalar or array_like of shape(M, )
        Latitude [radians].
    delta : scalar or array_like of shape(M, )
        Earth declination [radians].

    Returns
    -------
    ndarray
        Sunset hour angle [radians].

    """
    return np.arccos(-np.tan(lat) * np.tan(delta))


def _ra_daily(lat, doy, method='refet'):
    """Daily extraterrestrial radiation (Eq. 21)

    Parameters
    ----------
    lat : scalar or array_like of shape(M, )
        latitude [radians].
    doy : scalar or array_like of shape(M, )
        Day of year.
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.

    Returns
    -------
    ra : ndarray
        Daily extraterrestrial radiation [MJ m-2 d-1].

    Notes
    -----
    Equation in ASCE-EWRI 2005 uses a solar constant of ~1366.666... W m-2
    Equation in Duffie & Beckman (?) uses a solar constant of 1367 W m-2

    """
    delta = _delta(doy, method)
    omegas = _omega_sunset(lat, delta)
    theta = (omegas * np.sin(lat) * np.sin(delta) +
             np.cos(lat) * np.cos(delta) * np.sin(omegas))
    # print('{:>10s}: {:>8.3f}'.format('delta', float(delta)))
    # print('{:>10s}: {:>8.3f}'.format('omegas', float(omegas)))
    # print('{:>10s}: {:>8.3f}'.format('theta', float(theta)))
    # print('{:>10s}: {:>8.3f}'.format('dr', float(_dr(doy))))
    if method == 'asce':
        ra = (24. / math.pi) * 4.92 * _dr(doy) * theta
    else:
        ra = (24. / math.pi) * (1367 * 0.0036) * _dr(doy) * theta
    return ra


def _ra_hourly(lat, lon, doy, time_mid, method='refet'):
    """Hourly extraterrestrial radiation (Eq. 48)

    Parameters
    ----------
    lat : scalar or array_like of shape(M, )
        Latitude [radians].
    lon : scalar or array_like of shape(M, )
        Longitude [radians].
    doy : scalar or array_like of shape(M, )
        Day of year.
    time_mid : scalar or array_like of shape(M, )
        UTC time at midpoint of period [hours].
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.

    Returns
    -------
    ra : ndarray
        Hourly extraterrestrial radiation [MJ m-2 h-1].

    Notes
    -----
    Equation in ASCE-EWRI 2005 uses a solar constant of ~1366.666... W m-2
    Equation in Duffie & Beckman (?) uses a solar constant of 1367 W m-2
    """
    omega = _omega(_solar_time_rad(lon, time_mid, _seasonal_correction(doy)))
    delta = _delta(doy, method)
    omegas = _omega_sunset(lat, delta)

    # Solar time as start and end of period (Eqs. 53 & 54)
    # Modify omega1 and omega2 at sunrise and sunset (Eq. 56)
    omega1 = np.clip(omega - (math.pi / 24), -omegas, omegas)
    omega2 = np.clip(omega + (math.pi / 24), -omegas, omegas)
    omega1 = np.minimum(omega1, omega2)

    # Extraterrestrial radiation (Eq. 48)
    theta = (
        ((omega2 - omega1) * np.sin(lat) * np.sin(delta)) +
        (np.cos(lat) * np.cos(delta) * (np.sin(omega2) - np.sin(omega1))))
    if method == 'asce':
        ra = (12. / math.pi) * 4.92 * _dr(doy) * theta
    else:
        ra = (12. / math.pi) * (1367 * 0.0036) * _dr(doy) * theta
    return ra


def _rso_daily(ra, ea, pair, doy, lat):
    """Full daily clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ra : scalar or array_like of shape(M, )
        Extraterrestrial radiation [MJ m-2 d-1].
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].
    pair : scalar or array_like of shape(M, )
        Air pressure [kPa].
    doy : scalar or array_like of shape(M, )
        Day of year.
    lat : scalar or array_like of shape(M, )
        Latitude [rad].

    Returns
    -------
    rso : ndarray
        Daily clear sky solar radiation [MJ m-2 d-1]

    """
    # sin of the angle of the sun above the horizon (D.5 and Eq. 62)
    sin_beta_24 = np.sin(
        0.85 + 0.3 * lat * np.sin(_doy_fraction(doy) - 1.39) -
        0.42 * np.power(lat, 2))
    sin_beta_24 = np.maximum(sin_beta_24, 0.1)

    # Precipitable water
    w = _precipitable_water(pair, ea)

    # Clearness index for direct beam radiation (Eq. D.2)
    # Limit sin_beta >= 0.01 so that KB does not go undefined
    kb = (0.98 * np.exp((-0.00146 * pair) / sin_beta_24 -
                        0.075 * np.power((w / sin_beta_24), 0.4)))

    # Transmissivity index for diffuse radiation (Eq. D.4)
    kd = np.minimum(-0.36 * kb + 0.35, 0.82 * kb + 0.18)

    # print('{:>10s}: {:>8.3f}'.format('sin_beta_24', float(sin_beta_24)))
    # print('{:>10s}: {:>8.3f}'.format('w', float(w)))
    # print('{:>10s}: {:>8.3f}'.format('kb', float(kb)))
    # print('{:>10s}: {:>8.3f}'.format('kd', float(kd)))

    rso = ra * (kb + kd)
    return rso


def _rso_hourly(ra, ea, pair, doy, time_mid, lat, lon, method='refet'):
    """Full hourly clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ra : scalar or array_like of shape(M, )
        Extraterrestrial radiation [MJ m-2 h-1].
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].
    pair : scalar or array_like of shape(M, )
        Air pressure [kPa].
    doy : scalar or array_like of shape(M, )
        Day of year.
    time_mid : scalar or array_like of shape(M, )
        UTC time at midpoint of period [hours].
    lat : scalar or array_like of shape(M, )
        Latitude [rad].
    lon : scalar or array_like of shape(M, )
        Longitude [rad].
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        Passed through to declination calculation (_delta()).

    Returns
    -------
    rso : ndarray
        Hourly clear sky solar radiation [MJ m-2 h-1].

    """
    sc = _seasonal_correction(doy)
    omega = _omega(_solar_time_rad(lon, time_mid, sc))

    # sin of the angle of the sun above the horizon (D.6 and Eq. 62)
    delta = _delta(doy, method)
    sin_beta = (
        np.sin(lat) * np.sin(delta) +
        np.cos(lat) * np.cos(delta) * np.cos(omega))

    # Precipitable water
    w = _precipitable_water(pair, ea)

    # Clearness index for direct beam radiation (Eq. D.2)
    # Limit sin_beta >= 0.01 so that KB does not go undefined
    kt = 1.0
    kb = 0.98 * np.exp(
        (-0.00146 * pair) / (kt * np.maximum(sin_beta, 0.01)) -
        0.075 * np.power((w / np.maximum(sin_beta, 0.01)), 0.4))

    # Transmissivity index for diffuse radiation (Eq. D.4)
    kd = np.minimum(-0.36 * kb + 0.35, 0.82 * kb + 0.18)

    rso = ra * (kb + kd)
    return rso


def _rso_simple(ra, elev):
    """Simplified daily/hourly clear sky solar formulation (Eqs. 19 & 45)

    Parameters
    ----------
    ra : scalar or array_like of shape(M, )
        Extraterrestrial radiation [MJ m-2 d-1 or MJ m-2 h-1].
    elev : scalar or array_like of shape(M, )
        Elevation [m].

    Returns
    -------
    rso : ndarray
        Clear sky solar radiation [MJ m-2 d-1 or MJ m-2 h-1].

    """
    rso = (0.75 + 2E-5 * elev) * ra
    return rso


def _fcd_daily(rs, rso):
    """Daytime cloudiness fraction (Eq. 18)

    Parameters
    ----------
    rs : scalar or array_like of shape(M, )
        Measured solar radiation [MJ m-2 d-1].
    rso : scalar or array_like of shape(M, )
        Clear sky solar radiation [MJ m-2 d-1].

    Returns
    -------
    ndarray

    """
    return 1.35 * np.clip(rs / rso, 0.3, 1.0) - 0.35


def _fcd_hourly(rs, rso, doy, time_mid, lat, lon, method='refet'):
    """Cloudiness fraction (Eq. 45)

    Parameters
    ----------
    rs : array_like of shape(M, )
        Measured solar radiation [MJ m-2 h-1].
    rso : array_like of shape(M, )
        Clear sky solar radiation [MJ m-2 h-1].
    doy : scalar or array_like of shape(M, )
        Day of year.
    time_mid : scalar or array_like of shape(M, )
        UTC time at midpoint of period [hours].
    lat : scalar or array_like of shape(M, )
        Latitude [rad].
    lon : scalar or array_like of shape(M, )
        Longitude [rad].
    method : {'refet', 'asce'}, optional
        Calculation method:
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
        Passed through to declination calculation (_delta()).

    Returns
    -------
    ndarray

    """
    rs = np.array(rs, copy=True, ndmin=1).astype(np.float64)
    rso = np.array(rso, copy=True, ndmin=1).astype(np.float64)

    sc = _seasonal_correction(doy)
    delta = _delta(doy, method)
    omega = _omega(_solar_time_rad(lon, time_mid, sc))
    beta = np.arcsin(
        np.sin(lat) * np.sin(delta) +
        np.cos(lat) * np.cos(delta) * np.cos(omega))

    fcd = np.ones(rso.shape)
    fcd[rso > 0] = 1.35 * np.clip(rs[rso > 0] / rso[rso > 0], 0.3, 1) - 0.35

    # For now just set fcd to 1 for low sun angles
    # DEADBEEF - Still need to get daytime value of fcd when beta > 0.3
    # Get closest value in time (array space) when beta > 0.3
    fcd[beta < 0.3] = 1

    return fcd


def _rnl_daily(tmax, tmin, ea, fcd):
    """Daily net long-wave radiation  (Eq. 17)

    Parameters
    ----------
    tmax : scalar or array_like of shape(M, )
        Daily maximum air temperature [C].
    tmin : scalar or array_like of shape(M, )
        Daily minimum air temperature [C].
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].
    fcd : scalar or array_like of shape(M, )
        cloudiness fraction.

    Returns
    -------
    ndarray
        Daily net long-wave radiation [MJ m-2 d-1].

    """
    rnl = (
        4.901E-9 * fcd * (0.34 - 0.14 * np.sqrt(ea)) *
        0.5 * (np.power(tmax + 273.15, 4) + np.power(tmin + 273.15, 4)))
    return rnl


def _rnl_hourly(tmean, ea, fcd):
    """Hourly net long-wave radiation  (Eq. 44)

    Parameters
    ----------
    tmean : scalar or array_like of shape(M, )
        Mean hourly air temperature [C].
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].
    fcd : scalar or array_like of shape(M, )
        Cloudiness fraction.

    Returns
    -------
    ndarray
        Hourly net long-wave radiation [MJ m-2 h-1].

    """
    rnl = (
         2.042E-10 * fcd * (0.34 - 0.14 * np.sqrt(ea)) *
         np.power((tmean + 273.16), 4))
    return rnl


def _wind_height_adjust(uz, zw):
    """Wind speed at 2 m height based on full logarithmic profile (Eq. 33)

    Parameters
    ----------
    uz : scalar or array_like of shape(M, )
        Wind speed at measurement height [m/s].
    zw : scalar or array_like of shape(M, )
        Wind measurement height [m].

    Returns
    -------
    ndarray
        Wind speed at 2 m height [m/s].

    """
    return uz * 4.87 / np.log(67.8 * zw - 5.42)

import math

import numpy as np


def daily(tmin, tmax, ea, rs, uz, zw, elev, lat, doy,
          ref_type='etr', rso_type='full', rso=None):
    """ASCE Daily Standardized Reference Evapotranspiration (ET)

    Arguments
    ---------
    tmin : ndarray
        Minimum daily temperature [C].
    tmax : ndarray
        Maximum daily temperature [C].
    ea : ndarray
        Actual vapor pressure [kPa].
    rs : ndarray
        Incoming shortwave solar radiation [MJ m-2 day-1].
    uz : ndarray
        Wind speed [m/s].
    zw : float
        Wind speed height [m].
    elev : ndarray
        Elevation [m].
    lat : ndarray
        Latitude [radians].
    doy : ndarray
        Day of year.
    ref_type : {'eto', 'etr', 'grass', 'alfalfa', 'short', 'tall'}, optional
        Specifies which reference crop surface:

        * 'etr' -- Tall reference crop (default)
        * 'alfalfa' -- Tall reference crop
        * 'tall' -- Tall reference crop
        * 'eto' -- Short reference crop
        * 'grass' -- Short reference crop
        * 'short' -- Short reference crop

    rso_type : {'full' (default), 'simple', 'array'}, optional
        Specifies which clear sky solar radiation (Rso) model to use:

        * 'full' -- Full clear sky solar formulation
        * 'simple' -- Simplified clear sky solar formulation (Eq. 19)
        * 'array' -- Read Rso values from "rso" function parameter

    rso : float or None, optional
        Clear sky solar radiation [MJ m-2 day-1].
        Only needed if rso_type is 'array'.

    Returns
    -------
    etsz : ndarray
        Standardized reference ET [mm].

    Raises
    ------
    ValueError
        If "ref_type" or "rso_type" are invalid.

    Notes
    -----
    cn: 900 for ETo, 1600 for ETr
    cd: 0.34 for ETo, 0.38 for ETr
    Divide solar radiation values by 0.0864 to convert MJ m-2 day-1 to W m-2

    References
    ----------
    .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
        equation. ASCE-EWRI Standardization of Reference Evapotranspiration
        Task Committee Rep., ASCE Reston, Va.
        http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf
        http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf
    """

    # Convert all inputs to NumPy arrays
    tmin = np.array(tmin, copy=True, ndmin=1)
    tmax = np.array(tmax, copy=True, ndmin=1)
    ea = np.array(ea, copy=True, ndmin=1)
    rs = np.array(rs, copy=True, ndmin=1)
    uz = np.array(uz, copy=True, ndmin=1)
    elev = np.array(elev, copy=True, ndmin=1)
    lat = np.array(lat, copy=True, ndmin=1)

    # Check that latitudes are in radians
    if np.any(np.fabs(lat) > (0.5 * math.pi)):
        raise ValueError('latitudes must be in radians [-pi/2, pi/2]')

    if ref_type.lower() in ['eto', 'grass', 'short']:
        # Tall reference crop parameters
        cn, cd = 900, 0.34
    elif ref_type.lower() in ['etr', 'alfalfa', 'tall']:
        # Short reference crop parameters
        cn, cd = 1600, 0.38
    else:
        raise ValueError('ref_type must be "etr" or "eto"')

    # To match standardized form, psy is calculated from elevation based pair
    pair = _air_pressure(elev)
    psy = 0.000665 * pair

    # Vapor pressure
    tmean = 0.5 * (tmax + tmin)
    es_slope = (
        4098 * _sat_vapor_pressure(tmean) / (np.power((tmean + 237.3), 2)))
    es = 0.5 * (_sat_vapor_pressure(tmax) + _sat_vapor_pressure(tmin))

    # DEADBBEF - remove
    # Vapor pressure from RHmax and RHmin
    # ea = 0.5 * (es_tmin * rhmax + es_tmax * rhmin)

    # DEADBBEF - remove
    # Vapor pressure from specific humidity
    # To match standardized form, ea is calculated from elevation based pair
    # ea = _actual_vapor_pressure_func(q, pair)

    # Extraterrestrial radiation
    ra = _ra_daily(lat, doy)

    # Clear sky solar radiation
    if rso_type.lower() == 'full':
        # This is the full clear sky solar formulation
        rso = _rso_daily(ra, ea, pair, doy, lat)
    elif rso_type.lower() == 'simple':
        # Simplified clear sky solar formulation (Eq. 19)
        rso = _rso_simple(ra, elev)
    elif rso_type.lower() == 'array':
        # Use rso array passed to function
        pass
    else:
        raise ValueError('rso_type must be "simple", "full", or "array')

    # Cloudiness fraction (Eq. 18)
    fcd = 1.35 * np.clip(rs / rso, 0.3, 1.0) - 0.35

    # Net long-wave radiation
    rnl = _rnl_daily(tmax, tmin, ea, fcd)

    # Net radiation (Eqs. 15 and 16)
    rn = 0.77 * rs - rnl

    # Wind speed
    u2 = _wind_height_adjust(uz, zw)

    # Daily reference ET (Eq. 1)
    etsz = (
        (0.408 * es_slope * rn + (psy * cn * u2 * (es - ea) / (tmean + 273))) /
        (es_slope + psy * (cd * u2 + 1)))

    # print('\n{:>10s}: {:>8.3f}'.format('tmin', float(tmin)))
    # print('{:>10s}: {:>8.3f}'.format('tmax', float(tmax)))
    # print('{:>10s}: {:>8.3f}'.format('tmean', float(0.5 * (tmax + tmin))))
    # print('{:>10s}: {:>8.3f}'.format('rs', float(rs)))
    # print('{:>10s}: {:>8.3f}'.format('lat', float(lat)))
    # print('{:>10s}: {:>8.3f}'.format('pair', float(pair)))
    # print('{:>10s}: {:>8.3f}'.format('ea', float(ea)))
    # print('{:>10s}: {:>8.3f}'.format('es', float(es)))
    # print('{:>10s}: {:>8.3f}'.format('es_slope', float(es_slope)))
    # print('{:>10s}: {:>8.3f}'.format('ra', float(ra)))
    # print('{:>10s}: {:>8.3f}'.format('rs', float(rs)))
    # print('{:>10s}: {:>8.3f}'.format('rso', float(rso)))
    # print('{:>10s}: {:>8.3f}'.format('Fcd', float(fcd)))
    # print('{:>10s}: {:>8.3f}'.format('Rnl', float(rnl)))
    # print('{:>10s}: {:>8.3f}'.format('Rn', float(rn)))
    # print('{:>10s}: {:>8.3f}'.format('u2', float(u2)))
    # print('{:>10s}: {:>8.3f}'.format('etsz', float(etsz)))

    return etsz


def hourly(tmean, ea, rs, uz, zw, elev, lat, lon, doy, time, ref_type='etr'):
    """ASCE Hourly Standardized Reference Evapotranspiration (ET)

    .. warning:: Cloudiness fraction at night is not being computed correctly

    Arguments
    ---------
    tmean : array_like
        Average hourly temperature [C].
    ea : ndarray
        Actual vapor pressure [kPa].
    rs : array_like
        Shortwave solar radiation [MJ m-2 day-1].
    uz : array_like
        Wind speed [m/s].
    zw : float
        Wind speed measurement/estimated height [m].
    elev : array_like
        Elevation [m]
    lat : array_like
        Latitude [radians]
    lon : array_like
        Longitude [radians].
    doy : array_like
        Day of year.
    time : array_like
        UTC hour at start of time period.
    ref_type : {'eto', 'etr', 'grass', 'alfalfa', 'short', 'tall'}, optional
        Specifies which reference crop surface.

        * 'etr' -- Tall reference crop (default)
        * 'alfalfa' -- Tall reference crop
        * 'tall' -- Tall reference crop
        * 'eto' -- Short reference crop
        * 'grass' -- Short reference crop
        * 'short' -- Short reference crop

    Returns
    -------
    etsz : ndarray
        Standardized reference ET [mm].

    Raises
    ------
    ValueError
        If "ref_type" is invalid.

    Notes
    -----
    Divide solar radiation values by 0.0036 to convert MJ m-2 hr-1 to W m-2

    References
    ----------
    .. [1] ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
        equation. ASCE-EWRI Standardization of Reference Evapotranspiration
        Task Committee Rep., ASCE Reston, Va.
        http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf
        http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf
    """

    # Convert all inputs to NumPy arrays
    tmean = np.array(tmean, copy=True, ndmin=1)
    ea = np.array(ea, copy=True, ndmin=1)
    rs = np.array(rs, copy=True, ndmin=1)
    uz = np.array(uz, copy=True, ndmin=1)
    elev = np.array(elev, copy=True, ndmin=1)
    lat = np.array(lat, copy=True, ndmin=1)
    lon = np.array(lon, copy=True, ndmin=1)
    doy = np.array(doy, copy=True, ndmin=1)
    time = np.array(time, copy=True, ndmin=1)
    time_mid = time + 0.5

    # Check that latitude & longitude are in radians
    if np.any(np.fabs(lat) > (0.5 * math.pi)):
        raise ValueError('latitudes must be in radians [-pi/2, pi/2]')
    elif np.any(np.fabs(lon) > math.pi):
        raise ValueError('longitudes must be in radians [-pi, pi]')

    if ref_type.lower() in ['eto', 'grass', 'short']:
        # Short reference crop parameters
        cn_day = 37.0
        cd_day = 0.24
        g_rn_day = 0.1
        cn_night = 37.0
        cd_night = 0.96
        g_rn_night = 0.5
    elif ref_type.lower() in ['etr', 'alfalfa', 'tall']:
        # Tall reference crop parameters
        cn_day = 66.0
        cd_day = 0.25
        g_rn_day = 0.04
        cn_night = 66.0
        cd_night = 1.7
        g_rn_night = 0.2
    else:
        raise ValueError('ref_type must be "etr" or "eto"')

    # To match standardized form, psy is calculated from elevation based pair
    pair = _air_pressure(elev)
    psy = 0.000665 * pair
    es = _sat_vapor_pressure(tmean)
    es_slope = 4098 * es / np.power((tmean + 237.3), 2)

    # DEADBBEF - remove
    # Vapor pressure from specific humidity
    # To match standardized form, ea is calculated from elevation based pair
    # ea = _actual_vapor_pressure_func(q, pair)

    # Extraterrestrial radiation
    ra = _ra_hourly(lat, lon, doy, time_mid)

    # Simplified clear sky solar radiation
    # rso = _rso_simple(ra, elev)

    # Clear sky solar radiation
    rso = _rso_hourly(ra, ea, pair, doy, time_mid, lat, lon)

    # DEADBEEF - Move fcd calc to function
    # Cloudiness fraction (Eq. 45)
    # In IN2, "Beta" is computed for the start of the time period,
    #   but "SinBeta" is computed for the midpoint.
    # Beta (not SinBeta) is used for clamping fcd.
    sc = _seasonal_correction(doy)
    delta = _delta(doy)
    omega = _omega(_solar_time_rad(lon, time, sc))
    beta = np.arcsin(
        np.sin(lat) * np.sin(delta) +
        np.cos(lat) * np.cos(delta) * np.cos(omega))

    fcd = np.ones(beta.shape)
    fcd[rso > 0] = 1.35 * np.clip(rs[rso > 0] / rso[rso > 0], 0.3, 1) - 0.35

    # For now just set fcd to 1 for low sun angles
    # DEADBEEF - Still need to get daytime value of fcd when beta > 0.3
    # Get closest value in time (array space) when beta > 0.3
    fcd[beta < 0.3] = 1

    # Net long-wave radiation
    rnl = _rnl_hourly(tmean, ea, fcd)

    # Net radiation (Eqs. 42 and 43)
    rn = rs * 0.77 - rnl

    # Adjust coefficients for daytime/nighttime
    # Nighttime is defined as when Rn < 0 (pg 44)
    cn = np.zeros(rn.shape)
    cd = np.zeros(rn.shape)
    g_rn = np.zeros(rn.shape)
    cn[:] = cn_day
    cd[:] = cd_day
    g_rn[:] = g_rn_day
    rn_mask = rn < 0
    cn[rn_mask] = cn_night
    cd[rn_mask] = cd_night
    g_rn[rn_mask] = g_rn_night

    # Soil heat flux (Eqs. 65 and 66)
    g = rn * g_rn

    # Wind speed
    u2 = _wind_height_adjust(uz, zw)

    # Hourly reference ET (Eq. 1)
    etsz = (
        (0.408 * es_slope * (rn - g) + (psy * cn * u2 * (es - ea) / (tmean + 273))) /
        (es_slope + psy * (cd * u2 + 1)))

    # print('\n{:>10s}: {}'.format('tmean', float(tmean)))
    # print('{:>10s}: {}'.format('ea', float(ea)))
    # print('{:>10s}: {}'.format('rs', float(rs) / 0.0036))
    # print('{:>10s}: {}'.format('time (mid)', float(time_mid)))
    # print('{:>10s}: {}'.format('lat', float(lat)))
    # print('{:>10s}: {}'.format('lon', float(lon)))
    # print('{:>10s}: {}'.format('pair', float(pair)))
    # print('{:>10s}: {}'.format('es', float(es)))
    # print('{:>10s}: {}'.format('es_slope', float(es_slope)))
    # print('{:>10s}: {}'.format('ra', float(ra) / 0.0036))
    # print('{:>10s}: {}'.format('rs', float(rs) / 0.0036))
    # print('{:>10s}: {}'.format('rso', float(rso) / 0.0036))
    # print('{:>10s}: {}'.format('sc', float(sc)))
    # print('{:>10s}: {}'.format('delta', float(delta)))
    # print('{:>10s}: {}'.format('omega', float(omega)))
    # print('{:>10s}: {}'.format('beta', float(beta)))
    # print('{:>10s}: {}'.format('Fcd', float(fcd)))
    # print('{:>10s}: {}'.format('Rnl', float(rnl) / 0.0036))
    # print('{:>10s}: {}'.format('Rn', float(rn) / 0.0036))
    # print('{:>10s}: {}'.format('G', float(g) / 0.0036))
    # print('{:>10s}: {}'.format('u2', float(u2)))
    # print('{:>10s}: {}'.format('etsz', float(etsz)))

    return etsz


def _delta(doy):
    """Earth declination (Eq. 51)

    Parameters
    ----------
    doy : array_like
        Day of year.

    Returns
    -------
    ndarray
        Earth declination [radians].

    Notes
    -----
    Additional decimal places are needed to match RefET exactly.
    Eq 24 in report is: 0.409 * sin(doy_fraction_func(doy) - 1.39)

    """
    return 0.40928 * np.sin(_doy_fraction(doy) - 1.39435)


def _doy_fraction(doy):
    """Fraction of the DOY in the year (Eq. 50)

    Parameters
    ----------
    doy : array_like
        Day of year.

    Returns
    -------
    array_like
        DOY fraction [radians].

    """
    return doy * (2 * math.pi / 365.)


def _sat_vapor_pressure(temperature):
    """Saturation vapor pressure from temperature (Eq. 7)

    Parameters
    ----------
    temperature : array_like
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


def _actual_vapor_pressure(q, pair):
    """"Actual vapor pressure from specific humidity

    Parameters
    ----------
    q : array_like
        Specific humidity [kg/kg].
    pair : array_like
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
    ea : array_like
        Specific humidity [kPa].
    pair : array_like
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


def _air_pressure(elev, temperature=20):
    """Mean atmospheric pressure at station elevation (Eqs. 3 & 34)

    Parameters
    ----------
    elev : array_like
        Elevation [m].
    temperature : array_like, optional
        Temperature [C].

    Returns
    -------
    pair : ndarray
        Air pressure [kPa].

    Notes
    -----
    101.3 * (((293.15 - 0.0065 * elev) / 293.15) ** 5.26)

    In some older version of RefET, pair was calculated using the following
        101.3 * ((285 - 0.0065 * elev) / 285 ** 5.26)

    """
    pair = np.array(elev, copy=True, ndmin=1).astype(np.float64)
    pair *= -0.0065
    pair += (temperature + 273.15)
    pair /= (temperature + 273.15)
    np.power(pair, 5.26, out=pair)
    pair *= 101.3
    return pair


def _precipitable_water(pair, ea):
    """Precipitable water in the atmosphere (Eq. D.3)

    Parameters
    ----------
    pair : array_like
        Air pressure [kPa].
    ea : array_like
        Vapor pressure [kPa].

    Returns
    -------
    ndarray
        Precipitable water [mm].

    """
    return pair * 0.14 * ea + 2.1


def _dr(doy):
    """Inverse square of the Earth-Sun Distance (Eq. 50)

    Parameters
    ----------
    doy : array_like
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
    doy : ndarray
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
    lon : ndarray
        Longitude [radians].
    time_mid : ndarray
        UTC time at midpoint of period [hours].
    sc : array_like
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
    solar_time : ndarray
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
    lat : ndarray
        Latitude [radians].
    delta : ndarray
        Earth declination [radians].

    Returns
    -------
    ndarray
        Sunset hour angle [radians].

    """
    return np.arccos(-np.tan(lat) * np.tan(delta))


def _ra_daily(lat, doy):
    """Daily extraterrestrial radiation (Eq. 21)

    Parameters
    ----------
    lat : ndarray
        latitude [radians].
    doy : ndarray
        Day of year.

    Returns
    -------
    ra : ndarray
        Daily extraterrestrial radiation [MJ m-2 d-1].

    """
    delta = _delta(doy)
    omegas = _omega_sunset(lat, delta)
    theta = (omegas * np.sin(lat) * np.sin(delta) +
             np.cos(lat) * np.cos(delta) * np.sin(omegas))
    # print('{:>10s}: {:>8.3f}'.format('delta', float(delta)))
    # print('{:>10s}: {:>8.3f}'.format('omegas', float(omegas)))
    # print('{:>10s}: {:>8.3f}'.format('theta', float(theta)))
    # print('{:>10s}: {:>8.3f}'.format('dr', float(_dr(doy))))
    ra = (24. / math.pi) * 4.92 * _dr(doy) * theta
    return ra


def _ra_hourly(lat, lon, doy, time_mid):
    """Hourly extraterrestrial radiation (Eq. 48)

    Parameters
    ----------
    lat : array_like
        Latitude [radians].
    lon : array_like
        Longitude [radians].
    doy : array_like
        Day of year.
    time_mid : array_like
        UTC time at midpoint of period [hours].

    Returns
    -------
    ra : ndarray
        Hourly extraterrestrial radiation [MJ m-2 h-1].

    """
    omega = _omega(_solar_time_rad(lon, time_mid, _seasonal_correction(doy)))
    delta = _delta(doy)
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
    ra = (12. / math.pi) * 4.92 * _dr(doy) * theta
    return ra


def _rso_daily(ra, ea, pair, doy, lat):
    """Full daily clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ra : array_like
        Extraterrestrial radiation [MJ m-2 d-1].
    ea : array_like
        Actual vapor pressure [kPa].
    pair : array_like
        Air pressure [kPa].
    doy : array_like
        Day of year.
    lat : array_like
        Latitude [rad].

    Returns
    -------
    rso : ndarray
        Daily clear sky solar radiation [MJ m-2 d-1]

    """

    # sin of the angle of the sun above the horizon (D.5 and Eq. 62)
    sin_beta_24 = np.sin(
        0.85 + 0.3 * lat * np.sin(_doy_fraction(doy) - 1.39435) -
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


def _rso_hourly(ra, ea, pair, doy, time_mid, lat, lon):
    """Full hourly clear sky solar radiation formulation (Appendix D)

    Parameters
    ----------
    ra : array_like
        Extraterrestrial radiation [MJ m-2 h-1].
    ea : array_like
        Actual vapor pressure [kPa].
    pair : array_like
        Air pressure [kPa].
    doy : array_like
        Day of year.
    time_mid : array_like
        UTC time at midpoint of period [hours].
    lat : array_like
        Latitude [rad].
    lon : array_like
        Longitude [rad].

    Returns
    -------
    rso : ndarray
        Hourly clear sky solar radiation [MJ m-2 h-1].

    """
    sc = _seasonal_correction(doy)
    omega = _omega(_solar_time_rad(lon, time_mid, sc))

    # sin of the angle of the sun above the horizon (D.6 and Eq. 62)
    delta = _delta(doy)
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
    ra : array_like
        Extraterrestrial radiation [MJ m-2 d-1 or MJ m-2 h-1].
    elev : array_like
        Elevation [m].

    Returns
    -------
    rso : ndarray
        Clear sky solar radiation [MJ m-2 d-1 or MJ m-2 h-1].

    """
    rso = (0.75 + 2E-5 * elev) * ra
    return rso


def _rnl_daily(tmax, tmin, ea, fcd):
    """Daily net long-wave radiation  (Eq. 17)

    Parameters
    ----------
    tmax : array_like
        Daily maximum air temperature [C].
    tmin : array_like
        Daily minimum air temperature [C].
    ea : array_like
        Actual vapor pressure [kPa].
    fcd : array_like
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
    tmean : array_like):
        Mean hourly air temperature [C].
    ea : array_like):
        Actual vapor pressure [kPa].
    fcd : array_like):
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
    uz : array_like):
        Wind speed at measurement height [m/s].
    zw : array_like):
        Wind measurement height [m].

    Returns
    -------
    ndarray
        Wind speed at 2 m height [m/s].

    """
    return uz * 4.87 / np.log(67.8 * zw - 5.42)

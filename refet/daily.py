import math

import numpy as np

from . import calcs


def daily(tmin, tmax, ea, rs, uz, zw, elev, lat, doy, surface,
          method='refet', rso_type='full', rso=None):
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
    surface : {'eto', 'etr', 'grass', 'alfalfa', 'short', 'tall'}
        Specifies which reference crop surface.
        * 'etr', 'alfalfa', 'tall' -- Tall reference crop
        * 'eto', 'grass', 'short' -- Short reference crop
    method : {'refet', 'asce'}, optional
        Specifies which calculation method to use.
        * 'refet' -- Calculations will follow RefET software (default).
        * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
    rso_type : {'full' (default), 'simple', 'array'}, optional
        Specifies which clear sky solar radiation (Rso) model to use.
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
        If "surface" or "rso_type" are invalid.
        If latitude values are outside the range [-pi/2, pi/2].

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

    if method.lower() not in ['asce', 'refet']:
        raise ValueError('method must be "asce" or "refet"')

    if surface.lower() in ['eto', 'grass', 'short']:
        # Tall reference crop parameters
        cn, cd = 900, 0.34
    elif surface.lower() in ['etr', 'alfalfa', 'tall']:
        # Short reference crop parameters
        cn, cd = 1600, 0.38
    else:
        raise ValueError('surface must be "etr" or "eto"')

    # To match standardized form, psy is calculated from elevation based pair
    pair = calcs._air_pressure(elev, method)
    psy = 0.000665 * pair

    # Vapor pressure
    tmean = 0.5 * (tmax + tmin)
    es_slope = (
        4098 * calcs._sat_vapor_pressure(tmean) / (np.power((tmean + 237.3), 2)))
    es = 0.5 * (calcs._sat_vapor_pressure(tmax) + calcs._sat_vapor_pressure(tmin))

    # DEADBBEF - remove
    # Vapor pressure from RHmax and RHmin
    # ea = 0.5 * (es_tmin * rhmax + es_tmax * rhmin)

    # DEADBBEF - remove
    # Vapor pressure from specific humidity
    # To match standardized form, ea is calculated from elevation based pair
    # ea = _actual_vapor_pressure_func(q, pair)

    # Extraterrestrial radiation
    ra = calcs._ra_daily(lat, doy, method)

    # Clear sky solar radiation
    if rso_type.lower() == 'full':
        # This is the full clear sky solar formulation
        rso = calcs._rso_daily(ra, ea, pair, doy, lat)
    elif rso_type.lower() == 'simple':
        # Simplified clear sky solar formulation (Eq. 19)
        rso = calcs._rso_simple(ra, elev)
    elif rso_type.lower() == 'array':
        # Use rso array passed to function
        pass
    else:
        raise ValueError('rso_type must be "simple", "full", or "array')

    # Cloudiness fraction
    fcd = calcs._fcd_daily(rs, rso)

    # Net long-wave radiation
    rnl = calcs._rnl_daily(tmax, tmin, ea, fcd)

    # Net radiation (Eqs. 15 and 16)
    rn = 0.77 * rs - rnl

    # Wind speed
    u2 = calcs._wind_height_adjust(uz, zw)

    # Tmean units conversion factor
    # Check RefET to see if it is using 273 or 273.15
    if method.lower() == 'asce':
        c2k = 273.0
    elif method.lower() == 'refet':
        c2k = 273.0
        # c2k = 273.15

    # Daily reference ET (Eq. 1)
    etsz = (
        (0.408 * es_slope * rn + (psy * cn * u2 * (es - ea) / (tmean + c2k))) /
        (es_slope + psy * (cd * u2 + 1)))

    return etsz

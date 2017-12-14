import math

import numpy as np

from . import calcs


def hourly(tmean, ea, rs, uz, zw, elev, lat, lon, doy, time,
           ref_type='etr', asce_flag=False):
    """ASCE Hourly Standardized Reference Evapotranspiration (ET)

    .. warning:: Cloudiness fraction at night is not being computed correctly

    Arguments
    ---------
    tmean : array_like
        Average hourly temperature [C].
    ea : ndarray
        Actual vapor pressure [kPa].
    rs : array_like
        Shortwave solar radiation [MJ m-2 hr-1].
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
    asce_flag : bool
        if True, use  equations as defined in ASCE-EWRI 2005 [1].
        if False, use equations as defined in Ref-ET software.

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
    pair = calcs._air_pressure(elev, asce_flag)
    psy = 0.000665 * pair
    es = calcs._sat_vapor_pressure(tmean)
    es_slope = 4098 * es / np.power((tmean + 237.3), 2)

    # DEADBBEF - remove
    # Vapor pressure from specific humidity
    # To match standardized form, ea is calculated from elevation based pair
    # ea = _actual_vapor_pressure_func(q, pair)

    # Extraterrestrial radiation
    ra = calcs._ra_hourly(lat, lon, doy, time_mid, asce_flag)

    # Simplified clear sky solar radiation
    # rso = _rso_simple(ra, elev)

    # Clear sky solar radiation
    rso = calcs._rso_hourly(ra, ea, pair, doy, time_mid, lat, lon, asce_flag)

    # Cloudiness fraction
    # Intentionally not using time_mid to match Beta value in IN2 file
    # In IN2, "Beta" is computed for the start of the time period,
    #   but "SinBeta" is computed for the midpoint.
    # Beta (not SinBeta) is used for clamping fcd.
    fcd = calcs._fcd_hourly(rs, rso, doy, time, lat, lon, asce_flag)

    # Net long-wave radiation
    rnl = calcs._rnl_hourly(tmean, ea, fcd)

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
    u2 = calcs._wind_height_adjust(uz, zw)

    # Hourly reference ET (Eq. 1)
    etsz = (
        (0.408 * es_slope * (rn - g) + (psy * cn * u2 * (es - ea) / (tmean + 273))) /
        (es_slope + psy * (cd * u2 + 1)))

    # print('\n{:>10s}: {:>8.3f}'.format('tmean', float(tmean)))
    # print('{:>10s}: {:>8.3f}'.format('ea', float(ea)))
    # print('{:>10s}: {:>8.3f}'.format('rs', float(rs) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('time (mid)', float(time_mid)))
    # print('{:>10s}: {:>8.3f}'.format('lat', float(lat)))
    # print('{:>10s}: {:>8.3f}'.format('lon', float(lon)))
    # print('{:>10s}: {:>8.3f}'.format('pair', float(pair)))
    # print('{:>10s}: {:>8.3f}'.format('es', float(es)))
    # print('{:>10s}: {:>8.4f}'.format('es_slope', float(es_slope)))
    # print('{:>10s}: {:>8.3f}'.format('ra', float(ra) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('rs', float(rs) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('rso', float(rso) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('Fcd', float(fcd)))
    # print('{:>10s}: {:>8.3f}'.format('Rnl', float(rnl) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('Rn', float(rn) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('G', float(g) / 0.0036))
    # print('{:>10s}: {:>8.3f}'.format('u2', float(u2)))
    # print('{:>10s}: {:>8.3f}'.format('etsz', float(etsz)))

    return etsz

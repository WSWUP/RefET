import math

import numpy as np

from . import calcs


class Daily():
    def __init__(self, tmin, tmax, ea, rs, uz, zw, elev, lat, doy,
                 method='asce', rso_type=None, rso=None,
                 input_units={}, output_units={}):
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
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
            * 'refet' -- Calculations will follow RefET software.
        rso_type : {None (default), 'full' , 'simple', 'array'}, optional
            Specifies which clear sky solar radiation (Rso) model to use.
            * None -- Rso type will be determined from "method" parameter
            * 'full' -- Full clear sky solar formulation
            * 'simple' -- Simplified clear sky solar formulation
            * 'array' -- Read Rso values from "rso" function parameter
        rso : array_like or None, optional
            Clear sky solar radiation [MJ m-2 day-1] (the default is None).
            Only used if rso_type == 'array'.
        input_units : dict, optional
            Input unit types.
        output_units : dict, optional
            Output unit types.

        Returns
        -------
        etsz : ndarray
            Standardized reference ET [mm].

        Raises
        ------
        ValueError
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
        self.tmin = np.array(tmin, copy=True, ndmin=1)
        self.tmax = np.array(tmax, copy=True, ndmin=1)
        self.ea = np.array(ea, copy=True, ndmin=1)
        self.rs = np.array(rs, copy=True, ndmin=1)
        self.uz = np.array(uz, copy=True, ndmin=1)
        self.elev = np.array(elev, copy=True, ndmin=1)
        self.lat = np.array(lat, copy=True, ndmin=1)
        self.zw = zw
        self.doy = doy

        # Check that latitudes are in radians (after applying unit conversion)
        if np.any(np.fabs(lat) > (0.5 * math.pi)):
            raise ValueError('latitudes must be in radians [-pi/2, pi/2]')

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        if rso_type is None:
            pass
        elif rso_type.lower() not in ['simple', 'full', 'array']:
            raise ValueError('rso_type must be None, "simple", "full", or "array')
        elif rso_type.lower() in 'array':
            # Check that rso is an array
            pass

        # To match standardized form, pair is calculated from elevation
        self.pair = calcs._air_pressure(self.elev, method)

        self.psy = 0.000665 * self.pair

        self.tmean = 0.5 * (self.tmax + self.tmin)
        self.es_slope = calcs._es_slope(self.tmean, method)

        # Saturated vapor pressure
        self.es = 0.5 * (
            calcs._sat_vapor_pressure(self.tmax) +
            calcs._sat_vapor_pressure(self.tmin))

        # Vapor pressure deficit
        self.vpd = calcs._vpd(self.es, self.ea)

        # Extraterrestrial radiation
        self.ra = calcs._ra_daily(self.lat, self.doy, method)

        # Clear sky solar radiation
        # If rso_type is not set, use the method
        # If rso_type is set, use rso_type directly
        if rso_type is None :
            if method.lower() == 'asce':
                self.rso = calcs._rso_simple(self.ra, self.elev)
            elif method.lower() == 'refet':
                self.rso = calcs._rso_daily(
                    self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'simple':
            self.rso = calcs._rso_simple(self.ra, elev)
        elif rso_type.lower() == 'full':
            self.rso = calcs._rso_daily(
                self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'array':
            # Use rso array passed to function
            self.rso = rso

        # Cloudiness fraction
        self.fcd = calcs._fcd_daily(self.rs, self.rso)

        # Net long-wave radiation
        self.rnl = calcs._rnl_daily(self.tmax, self.tmin, self.ea, self.fcd)

        # Net radiation (Eqs. 15 and 16)
        self.rn = 0.77 * self.rs - self.rnl

        # Wind speed
        self.u2 = calcs._wind_height_adjust(self.uz, self.zw)

    def eto(self):
        """Grass reference surface"""
        self.cn = 900
        self.cd = 0.34
        return self._etsz()

    def etr(self):
        """Alfalfa reference surface"""
        self.cn = 1600
        self.cd = 0.38
        return self._etsz()

    def _etsz(self):
        """Daily reference ET (Eq. 1)"""
        return (
            (0.408 * self.es_slope * self.rn +
             (self.psy * self.cn * self.u2 * self.vpd / (self.tmean + 273))) /
            (self.es_slope + self.psy * (self.cd * self.u2 + 1)))

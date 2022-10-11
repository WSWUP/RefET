import math

import numpy as np

from . import calcs
from . import units


class Daily():
    def __init__(self, tmin, tmax, ea, rs, uz, zw, elev, lat, doy,
                 method='asce', rso_type=None, rso=None, input_units={}):
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
            Wind speed [m s-1].
        zw : float
            Wind speed height [m].
        elev : ndarray
            Elevation [m].
        lat : ndarray
            Latitude [degrees].
        doy : ndarray
            Day of year.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1]_ equations.
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

        Returns
        -------
        etsz : ndarray
            Standardized reference ET [mm].

        Notes
        -----
        cn: 900 for ETo, 1600 for ETr
        cd: 0.34 for ETo, 0.38 for ETr

        The Langleys to MJ m-2 conversion factor is the value used in the RefET
        program, although there are other factors that could be applied:
        https://www.aps.org/policy/reports/popa-reports/energy/units.cfm

        References
        ----------
        .. [1] ASCE-EWRI (2005). The ASCE standardized reference
           evapotranspiration equation. ASCE-EWRI Standardization of Reference
           Evapotranspiration Task Committee Rep., ASCE Reston, Va.

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

        # Unit conversion
        for v, unit in input_units.items():
            setattr(
                self, v,
                units.convert(getattr(self, v), v, unit, timestep='daily')
            )

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        if rso_type is None:
            pass
        elif rso_type.lower() not in ['simple', 'full', 'array']:
            raise ValueError('rso_type must be None, "simple", "full", or "array')
        elif rso_type.lower() in 'array':
            # Check that rso is an array
            pass

        # The input angles are converted to degrees by default in units.convert
        # They need to be converted back to radians for the calc functions
        # This is a little roundabout but was done to since the user is most
        #   likely using latitude values that are in degrees and would not be
        #   expecting the default units to be radians
        self.lat *= (math.pi / 180.0)

        # To match standardized form, pair is calculated from elevation
        self.pair = calcs._air_pressure(self.elev, method)

        # Psychrometric constant (Eq. 4)
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

        # Net radiation
        self.rn = calcs._rn_daily(self.rs, self.rnl)

        # Soil heat flux
        self.g = 0

        # Wind speed
        self.u2 = calcs._wind_height_adjust(self.uz, self.zw)

    def etsz(self, surface):
        """Standardized reference ET

        Parameters
        ----------
        surface : {'alfalfa', 'etr', 'tall', 'grass', 'eto', 'short'}
            Reference surface type.

        Returns
        -------
        ndarray

        """
        if surface.lower() in ['alfalfa', 'etr', 'tall']:
            return self.etr()
        elif surface.lower() in ['grass', 'eto', 'short']:
            return self.eto()
        else:
            raise ValueError(f'unsupported surface type: {surface}')

    def eto(self):
        """Grass reference surface"""
        self.cn = 900
        self.cd = 0.34
        return calcs._etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd)

    def etr(self):
        """Alfalfa reference surface"""
        self.cn = 1600
        self.cd = 0.38
        return calcs._etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd)

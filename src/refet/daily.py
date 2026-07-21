import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import calcs
from . import units


class Daily():
    def __init__(
            self,
            tmin: ArrayLike,
            tmax: ArrayLike,
            rs: ArrayLike,
            uz: ArrayLike,
            zw: float,
            elev: ArrayLike,
            lat: ArrayLike,
            doy: ArrayLike,
            ea: ArrayLike | None = None,
            tdew: ArrayLike | None = None,
            method: str = 'asce',
            rso_type: str | None = None,
            rso: ArrayLike | None = None,
            input_units: dict[str, str] | None = None,
            q: ArrayLike | None = None,
            rh_min: ArrayLike | None = None,
            rh_max: ArrayLike | None = None,
        ):
        """ASCE Daily Standardized Reference Evapotranspiration (ET)

        Arguments
        ---------
        tmin : ndarray
            Minimum daily temperature [C].
        tmax : ndarray
            Maximum daily temperature [C].
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
        ea : ndarray, optional
            Actual vapor pressure [kPa].  One of ea, tdew, q, or the rh_min/
            rh_max pair must be set; when several are given the first in
            that order is used (the ASCE-EWRI 2005 preference hierarchy).
        tdew : ndarray, optional
            Mean daily dew point temperature [C].
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 equations.
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
        q : ndarray, optional
            Specific humidity [kg kg-1].
        rh_min : ndarray, optional
            Minimum daily relative humidity [percent].
        rh_max : ndarray, optional
            Maximum daily relative humidity [percent].

        Notes
        -----
        Call eto(), etr(), or etsz(surface) to compute the standardized
        reference ET [mm], and results() for the intermediate terms.

        cn: 900 for ETo, 1600 for ETr
        cd: 0.34 for ETo, 0.38 for ETr

        The Langleys to MJ m-2 conversion factor is the value used in the RefET
        program, although there are other factors that could be applied:
        https://www.aps.org/policy/reports/popa-reports/energy/units.cfm

        References
        ----------
        ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration
        equation. ASCE-EWRI Standardization of Reference Evapotranspiration
        Task Committee Rep., ASCE Reston, Va.

        """
        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        # Convert all inputs to NumPy arrays
        self.tmin = np.array(tmin, copy=True, ndmin=1)
        self.tmax = np.array(tmax, copy=True, ndmin=1)
        self.rs = np.array(rs, copy=True, ndmin=1)
        self.uz = np.array(uz, copy=True, ndmin=1)
        self.elev = np.array(elev, copy=True, ndmin=1)
        self.lat = np.array(lat, copy=True, ndmin=1)
        self.zw = zw
        self.doy = np.array(doy, copy=True, ndmin=1)

        # Vapor pressure inputs, in the ASCE-EWRI 2005 preference order:
        # measured ea, then tdew, then specific humidity, then min/max RH
        if (rh_min is None) != (rh_max is None):
            raise ValueError('rh_min and rh_max must be set together')
        self.ea = np.array(ea, copy=True, ndmin=1) if ea is not None else None
        self.tdew = np.array(tdew, copy=True, ndmin=1) if tdew is not None else None
        self.q = np.array(q, copy=True, ndmin=1) if q is not None else None
        self.rh_min = np.array(rh_min, copy=True, ndmin=1) if rh_min is not None else None
        self.rh_max = np.array(rh_max, copy=True, ndmin=1) if rh_max is not None else None
        if all(v is None for v in [ea, tdew, q, rh_min]):
            raise ValueError(
                'One of "ea", "tdew", "q", or the "rh_min"/"rh_max" pair '
                'must be set')

        # Unit conversions
        if input_units is None:
            input_units = {}
        for v, unit in input_units.items():
            setattr(self, v, units.convert(getattr(self, v), v, unit, timestep='daily'))

        # Compute Ea after handling unit conversions so that the humidity
        # inputs are in their default units
        if self.ea is None:
            if self.tdew is not None:
                self.ea = calcs.sat_vapor_pressure(self.tdew)
            elif self.q is not None:
                self.ea = calcs.actual_vapor_pressure(
                    self.q, calcs.air_pressure(self.elev, method))
            else:
                self.ea = calcs.ea_from_rh_daily(
                    self.tmin, self.tmax, self.rh_min, self.rh_max)

        # Rso
        if rso_type is not None:
            if rso_type.lower() not in ['simple', 'full', 'array']:
                raise ValueError('rso_type must be None, "simple", "full", or "array"')
            if rso_type.lower() == 'array' and rso is None:
                raise ValueError('rso must be set when rso_type is "array"')

        # The input angles are converted to degrees by default in units.convert
        # They need to be converted back to radians for the calc functions
        # This is a little roundabout but was done to since the user is most
        #   likely using latitude values that are in degrees and would not be
        #   expecting the default units to be radians
        self.lat *= (math.pi / 180.0)

        # Mean daily air temperature
        self.tmean = 0.5 * (self.tmax + self.tmin)

        # To match standardized form, pair is calculated from elevation
        self.pair = calcs.air_pressure(self.elev, method)

        # Psychrometric constant (Eq. 4)
        self.psy = 0.000665 * self.pair

        # Slope of the saturation vapor pressure-temperature curve
        self.es_slope = calcs.es_slope(self.tmean, method)

        # Saturation vapor pressure
        self.es = 0.5 * (calcs.sat_vapor_pressure(self.tmax) + calcs.sat_vapor_pressure(self.tmin))

        # Vapor pressure deficit
        self.vpd = calcs.vpd(self.es, self.ea)

        # Extraterrestrial radiation
        self.ra = calcs.ra_daily(self.lat, self.doy, method)

        # Clear sky solar radiation
        # If rso_type is not set, use the method
        # If rso_type is set, use rso_type directly
        if rso_type is None:
            if method.lower() == 'asce':
                self.rso = calcs.rso_simple(self.ra, self.elev)
            elif method.lower() == 'refet':
                self.rso = calcs.rso_daily(self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'simple':
            self.rso = calcs.rso_simple(self.ra, self.elev)
        elif rso_type.lower() == 'full':
            self.rso = calcs.rso_daily(self.ra, self.ea, self.pair, self.doy, self.lat)
        elif rso_type.lower() == 'array':
            # Use rso array passed to function
            self.rso = np.array(rso, copy=True, ndmin=1)

        # Cloudiness fraction
        self.fcd = calcs.fcd_daily(self.rs, self.rso)

        # Net long-wave radiation
        self.rnl = calcs.rnl_daily(self.tmax, self.tmin, self.ea, self.fcd)

        # Net radiation
        self.rn = calcs.rn_daily(self.rs, self.rnl)

        # Soil heat flux
        self.g = 0

        # Wind speed
        self.u2 = calcs.wind_height_adjust(self.uz, self.zw)

    def results(self) -> dict[str, NDArray[np.float64]]:
        """Input and intermediate terms of the daily calculation

        Returns
        -------
        dict of ndarray
            Terms computed by the constructor, keyed by name: tmin, tmax,
            tmean, tdew (only when given), ea, es, es_slope, vpd, rs, pair,
            psy, ra, rso, fcd, rnl, rn, u2.

        Notes
        -----
        The surface dependent terms (cn, cd) are set by the eto() and etr()
        methods and are not included.

        """
        keys = [
            'tmin', 'tmax', 'tmean', 'tdew', 'ea', 'es', 'es_slope', 'vpd',
            'rs', 'pair', 'psy', 'ra', 'rso', 'fcd', 'rnl', 'rn', 'u2',
        ]
        return {k: getattr(self, k) for k in keys if getattr(self, k) is not None}

    def etsz(self, surface: str) -> NDArray[np.float64]:
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

    def eto(self) -> NDArray[np.float64]:
        """Grass reference surface"""
        self.cn = 900
        self.cd = 0.34
        return calcs.etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd
        )

    def etr(self) -> NDArray[np.float64]:
        """Alfalfa reference surface"""
        self.cn = 1600
        self.cd = 0.38
        return calcs.etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd
        )

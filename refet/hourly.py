import math

import numpy as np

from . import calcs


class Hourly():
    def __init__(self, tmean, ea, rs, uz, zw, elev, lat, lon, doy, time,
                 method='asce', input_units={}):
        """ASCE Hourly Standardized Reference Evapotranspiration (ET)

        Arguments
        ---------
        tmean : ndarray
            Average hourly temperature [C].
        ea : ndarray
            Actual vapor pressure [kPa].
        rs : ndarray
            Shortwave solar radiation [MJ m-2 hr-1].
        uz : ndarray
            Wind speed [m/s].
        zw : float
            Wind speed measurement/estimated height [m].
        elev : ndarray
            Elevation [m]
        lat : ndarray
            Latitude [radians]
        lon : ndarray
            Longitude [radians].
        doy : ndarray
            Day of year.
        time : ndarray
            UTC hour at start of time period.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1] equations.
            * 'refet' -- Calculations will follow RefET software.
        input_units : dict, optional
            Input unit types.

        Returns
        -------
        etsz : ndarray
            Standardized reference ET [mm].

        Raises
        ------
        ValueError
            If latitude values are outside the range [-pi/2, pi/2].
            If longitude values are outside the range [-pi, pi].

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
        self.tmean = np.array(tmean, copy=True, ndmin=1)
        self.ea = np.array(ea, copy=True, ndmin=1)
        self.rs = np.array(rs, copy=True, ndmin=1)
        self.uz = np.array(uz, copy=True, ndmin=1)
        self.elev = np.array(elev, copy=True, ndmin=1)
        self.lat = np.array(lat, copy=True, ndmin=1)
        self.lon = np.array(lon, copy=True, ndmin=1)
        self.doy = np.array(doy, copy=True, ndmin=1)
        self.time = np.array(time, copy=True, ndmin=1)
        self.time_mid = self.time + 0.5
        self.zw = zw
        self.doy = doy

        # Unit conversion
        for variable, unit in input_units.items():
            # Check input unit types
            if unit == '':
                continue
            elif unit.lower() in [
                    'c', 'celsius', 'm', 'meter', 'meters',
                    'rad', 'radians', 'kpa', 'm s-1', 'm/s',
                    'mj m-2 day-1', 'mj m-2 d-1']:
                continue
            elif unit.lower() not in [
                    'k', 'kelvin', 'f', 'fahrenheit',
                    'deg', 'degree', 'degrees',
                    'ft', 'feet', 'mph']:
                raise ValueError('unsupported unit conversion for {} {}'.format(
                    variable, unit))

            # Convert input values to expected units
            if variable == 'tmean':
                if unit.lower() in ['f', 'fahrenheit']:
                    self.tmean -= 32
                    self.tmean *= (5.0 / 9)
                elif unit.lower() in ['k', 'kelvin']:
                    self.tmean -= 273.15
            elif variable == 'uz':
                if unit.lower() in ['mph']:
                    self.uz *= 0.44704
            elif variable == 'zw':
                if unit.lower() in ['ft', 'feet']:
                    self.zw *= 0.3048
            elif variable == 'elev':
                if unit.lower() in ['ft', 'feet']:
                    self.elev *= 0.3048
            elif variable == 'lat':
                if unit.lower() in ['deg', 'degree', 'degrees']:
                    self.lat *= (math.pi / 180)
            elif variable == 'lon':
                if unit.lower() in ['deg', 'degree', 'degrees']:
                    self.lon *= (math.pi / 180)

        # Check that latitude & longitude are in radians
        if np.any(np.fabs(self.lat) > (0.5 * math.pi)):
            raise ValueError('latitudes must be in radians [-pi/2, pi/2]')
        elif np.any(np.fabs(self.lon) > math.pi):
            raise ValueError('longitudes must be in radians [-pi, pi]')

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        # To match standardized form, psy is calculated from elevation based pair
        self.pair = calcs._air_pressure(self.elev, method)
        self.psy = 0.000665 * self.pair
        self.es = calcs._sat_vapor_pressure(self.tmean)
        self.es_slope = calcs._es_slope(self.tmean, method)

        # Vapor pressure deficit
        self.vpd = self.es - self.ea
        # self.vpd = calcs._vpd(self.es, ea)

        # Extraterrestrial radiation
        self.ra = calcs._ra_hourly(
            self.lat, self.lon, self.doy, self.time_mid, method)

        # Clear sky solar radiation
        if method == 'asce':
            self.rso = calcs._rso_simple(self.ra, self.elev)
        elif method == 'refet':
            self.rso = calcs._rso_hourly(
                self.ra, self.ea, self.pair, self.doy, self.time_mid, self.lat,
                self.lon, method)

        # Cloudiness fraction
        # Intentionally not using time_mid to match Beta value in IN2 file
        # In IN2, "Beta" is computed for the start of the time period,
        #   but "SinBeta" is computed for the midpoint.
        # Beta (not SinBeta) is used for clamping fcd.
        self.fcd = calcs._fcd_hourly(
            self.rs, self.rso, self.doy, self.time, self.lat, self.lon, method)

        # Net long-wave radiation
        self.rnl = calcs._rnl_hourly(self.tmean, self.ea, self.fcd)

        # Net radiation (Eqs. 42 and 43)
        self.rn =  0.77 * self.rs - self.rnl

        # Soil heat flux (Eqs. 65 and 66)
        # self.g = self.rn * g_rn

        # Wind speed
        self.u2 = calcs._wind_height_adjust(self.uz, self.zw)

    def _etsz(self):
        """Hourly reference ET (Eq. 1)"""
        return (
            (0.408 * self.es_slope * (self.rn - self.g) +
             (self.psy * self.cn * self.u2 * self.vpd / (self.tmean + 273))) /
            (self.es_slope + self.psy * (self.cd * self.u2 + 1)))

    def eto(self):
        """Short (grass) reference surface"""
        self.cn = 37.0

        # Adjust coefficients for daytime/nighttime
        # Nighttime is defined as when Rn < 0 (pg 44)
        self.cd = np.full(self.rn.shape, 0.24)
        self.g_rn = np.full(self.rn.shape, 0.1)
        self.cd[self.rn < 0] = 0.96
        self.g_rn[self.rn < 0] = 0.5

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn * self.g_rn

        return self._etsz()

    def etr(self):
        """Tall (alfalfa) reference surface"""
        self.cn = 66.0

        # Adjust coefficients for daytime/nighttime
        # Nighttime is defined as when Rn < 0 (pg 44)
        self.cd = np.full(self.rn.shape, 0.25)
        self.g_rn = np.full(self.rn.shape, 0.04)
        self.cd[self.rn < 0] = 1.7
        self.g_rn[self.rn < 0] = 0.2

        # Soil heat flux (Eqs. 65 and 66)
        self.g = self.rn * self.g_rn

        return self._etsz()
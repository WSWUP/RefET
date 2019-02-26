import math

import numpy as np

from . import calcs


class Hourly():
    def __init__(self, tmean, ea, rs, uz, zw, elev, lat, lon, doy, time,
                 method='asce', input_units={}):
        """ASCE Hourly Standardized Reference Evapotranspiration (ET)

        .. warning:: Cloudiness fraction at night is not being computed per [1]_

        Arguments
        ---------
        tmean : ndarray
            Average hourly temperature [C].
        ea : ndarray
            Actual vapor pressure [kPa].
        rs : ndarray
            Shortwave solar radiation [MJ m-2 hr-1].
        uz : ndarray
            Wind speed [m s-1].
        zw : float
            Wind speed measurement/estimated height [m].
        elev : ndarray
            Elevation [m]
        lat : ndarray
            Latitude [degrees]
        lon : ndarray
            Longitude [degrees].
        doy : ndarray
            Day of year.
        time : ndarray
            UTC hour at start of time period.
        method : {'asce' (default), 'refet'}, optional
            Specifies which calculation method to use.
            * 'asce' -- Calculations will follow ASCE-EWRI 2005 [1]_ equations.
            * 'refet' -- Calculations will follow RefET software.
        input_units : dict, optional
            Input unit types.

        Returns
        -------
        etsz : ndarray
            Standardized reference ET [mm].

        Notes
        -----
        The Langleys to MJ m-2 conversion factor is the value used in the RefET
        program, although there are other factors that could be applied:
        https://www.aps.org/policy/reports/popa-reports/energy/units.cfm

        References
        ----------
        .. [1] ASCE-EWRI (2005). The ASCE standardized reference
           evapotranspiration equation. ASCE-EWRI Standardization of Reference
           Evapotranspiration Task Committee Rep., ASCE Reston, Va.
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
                    'c', 'celsius',
                    'mj m-2 hour-1', 'mj m-2 h-1',
                    'kpa',
                    'm s-1', 'm/s',
                    'm', 'meter', 'meters',
                    'deg', 'degree', 'degrees']:
                continue
            elif unit.lower() not in [
                    'k', 'kelvin', 'f', 'fahrenheit',
                    'pa',
                    'langleys', 'w m-2', 'w/m2',
                    'mph',
                    'ft', 'feet',
                    'rad', 'radian', 'radians']:
                raise ValueError('unsupported unit conversion for {} {}'.format(
                    variable, unit))

            # Convert input values to expected units
            if variable == 'tmean':
                if unit.lower() in ['f', 'fahrenheit']:
                    self.tmean -= 32
                    self.tmean *= (5.0 / 9)
                elif unit.lower() in ['k', 'kelvin']:
                    self.tmean -= 273.15
            elif variable == 'ea':
                if unit.lower() in ['pa']:
                    self.ea /= 1000.0
            elif variable == 'rs':
                if unit.lower() in ['langleys']:
                    self.rs *= 0.041868
                elif unit.lower() in ['w m-2', 'w/m2']:
                    self.rs *= 0.0036
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
                if unit.lower() in ['rad', 'radian', 'radians']:
                    # This is a little backwards but convert to degrees so that
                    # it can be converted to radians below.  This is done so
                    # that not setting the value will default to degrees.
                    self.lat *= (180.0 / math.pi)
            elif variable == 'lon':
                if unit.lower() in ['rad', 'radian', 'radians']:
                    # This is a little backwards but convert to degrees so that
                    # it can be converted to radians below.  This is done so
                    # that not setting the value will default to degrees.
                    self.lon *= (180.0 / math.pi)

        if method.lower() not in ['asce', 'refet']:
            raise ValueError('method must be "asce" or "refet"')

        # Convert latitude/longitude since calcs functions are expecting radians
        self.lat *= (math.pi / 180.0)
        self.lon *= (math.pi / 180.0)

        # To match standardized form, psy is calculated from elevation based pair
        self.pair = calcs._air_pressure(self.elev, method)

        # Psychrometric constant (Eq. 35)
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

        # Net radiation
        self.rn = calcs._rn_hourly(self.rs, self.rnl)

        # Soil heat flux (Eqs. 65 and 66)
        # self.g = self.rn * g_rn

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
            raise ValueError('unsupported surface type: {}'.format(surface))

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

        return calcs._etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd)

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

        return calcs._etsz(
            rn=self.rn, g=self.g, tmean=self.tmean, u2=self.u2, vpd=self.vpd,
            es_slope=self.es_slope, psy=self.psy, cn=self.cn, cd=self.cd)

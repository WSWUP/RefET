import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


DEFAULT_UNITS = [
    'c', 'celsius',
    'mj m-2 d-1', 'mj m-2 day-1',
    'mj m-2 h-1', 'mj m-2 hr-1', 'mj m-2 hour-1',
    'kpa',
    'm s-1', 'm/s',
    'm', 'meter', 'meters',
    'deg', 'degree', 'degrees',
]
RS_HOURLY_UNITS = [
    'w m-2 h-1', 'w m-2 hr-1', 'w m-2 hour-1', 'w/m2/hr', 'w/m2/hour'
]
RS_DAILY_UNITS = ['w m-2 d-1', 'w m-2 day-1', 'w/m2/d', 'w/m2/day']
SUPPORTED_UNITS = [
    'k', 'kelvin',
    'f', 'fahrenheit',
    'pa',
    'langleys',
    'w m-2', 'w/m2',
    'mph',
    'ft', 'feet',
    'rad', 'radian', 'radians',
]
SUPPORTED_UNITS.extend(RS_HOURLY_UNITS)
SUPPORTED_UNITS.extend(RS_DAILY_UNITS)


def convert(
        values: ArrayLike, variable: str | None, unit: str,
        timestep: str | None = None) -> FloatArray | float:
    """Unit conversion function

    Args:
        values: ndarray
        variable: str
        unit: str
        timestep: str, optional
            Must be set when converting Rs W/m2 values.
            Choices are "daily" or "hourly".

    Returns:
        ndarray

    Notes:
        The input is never modified in place; converted values are returned
        as a new object. Integer inputs are converted to floats as needed.

    """
    if unit == '':
        return values
    elif unit.lower() in DEFAULT_UNITS:
        return values
    elif unit.lower() not in SUPPORTED_UNITS:
        raise ValueError(f'unsupported unit conversion for {variable} {unit}')

    # Convert input values to expected units
    # Plain (non in-place) arithmetic so that the caller's array is not
    #   mutated and integer inputs are safely promoted to floats
    if variable in ['tmean', 'tmin', 'tmax', 'tdew']:
        if unit.lower() in ['f', 'fahrenheit']:
            values = (values - 32) * (5.0 / 9)
        elif unit.lower() in ['k', 'kelvin']:
            values = values - 273.15
    elif variable == 'ea':
        if unit.lower() in ['pa']:
            values = values / 1000.0
    elif variable == 'rs':
        if unit.lower() in ['langleys']:
            values = values * 0.041868
        elif unit.lower() in ['w m-2', 'w/m2']:
            if timestep is None or timestep.lower() not in ['daily', 'hourly']:
                raise ValueError(f'unsupported rs timestep parameter: {timestep}')
            elif timestep.lower() == 'daily':
                values = values * 0.0864
            else:
                values = values * 0.0036
        elif unit.lower() in RS_DAILY_UNITS:
            values = values * 0.0864
        elif unit.lower() in RS_HOURLY_UNITS:
            values = values * 0.0036
    elif variable == 'uz':
        if unit.lower() in ['mph']:
            values = values * 0.44704
    elif variable in ['zw', 'elev']:
        if unit.lower() in ['ft', 'feet']:
            values = values * 0.3048
    elif variable in ['lat', 'lon']:
        if unit.lower() in ['rad', 'radian', 'radians']:
            # This is a little backwards but convert to degrees so that
            # it can be converted to radians below.  This is done so
            # that not setting the value will default to degrees.
            values = values * (180.0 / math.pi)

    return values


def deg2rad(deg: ArrayLike) -> FloatArray | float:
    """Convert degrees to radians"""
    return deg * math.pi / 180.0


def rad2deg(rad: ArrayLike) -> FloatArray | float:
    """Convert radians to degrees"""
    return rad * 180.0 / math.pi


def c2f(c: ArrayLike) -> FloatArray | float:
    """Convert Celsius to Fahrenheit"""
    return c * (9.0 / 5) + 32


def f2c(f: ArrayLike) -> FloatArray | float:
    """Convert Fahrenheit to Celsius"""
    return (f - 32) * (5.0 / 9)

import math

import pytest

import src.refet.units as units


def test_deg2rad(d=30, r=(math.pi / 6)):
    assert units._deg2rad(d) == pytest.approx(r)


def test_rad2deg(r=(math.pi / 6), d=30):
    assert units._rad2deg(r) == pytest.approx(d)


def test_c2f(c=20, f=68):
    assert units._c2f(c) == f


def test_f2c(f=68, c=20):
    assert units._f2c(f) == c


def test_convert_positional_inputs():
    # Check that the input parameter names and order don't change
    assert (units.convert(212, 'tmean', 'F') ==
            units.convert(values=212, variable='tmean', unit='F'))
    assert (units.convert(1, 'Rs', 'W/m2', 'daily') ==
            units.convert(values=1, variable='Rs', unit='W/m2', timestep='daily'))


@pytest.mark.parametrize('unit', units.DEFAULT_UNITS)
def test_convert_no_conversion(unit):
    assert 2 == units.convert(values=2, variable=None, unit=unit)


@pytest.mark.parametrize('unit', ['farhenheit'])
def test_convert_unsupported_unit(unit):
    with pytest.raises(ValueError):
        units.convert(values=2, variable=None, unit=unit)


@pytest.mark.parametrize('timestep', ['monthly'])
def test_convert_unsupported_timestep(timestep):
    with pytest.raises(ValueError):
        units.convert(values=2, variable='rs', unit='w/m2', timestep=timestep)


@pytest.mark.parametrize(
    'variable, unit, value, expected',
    [
        ['tmean', 'F', 32, 0],
        ['tmax', 'Fahrenheit', 212, 100],
        ['tmean', 'K', 273.15, 0],
        ['tmax', 'Kelvin', 373.15, 100],
    ]
)
def test_convert_temperature(variable, unit, value, expected):
    assert expected == units.convert(value, variable, unit)


@pytest.mark.parametrize(
    'variable, unit, value, expected',
    [
        ['ea', 'Pa', 1000, 1],
    ]
)
def test_convert_pressure(variable, unit, value, expected):
    assert expected == units.convert(value, variable, unit)


@pytest.mark.parametrize(
    'variable, unit, timestep, value, expected',
    [
        ['rs', 'langleys', None, 1, 0.041868],
        ['rs', 'W/m2', 'daily', 1, 0.0864],
        ['rs', 'W m-2', 'hourly', 1, 0.0036],
        # Assuming we only need to test one of the supported Rs unit combinations
        ['rs', 'w m-2 d-1', None, 1, 0.0864],
        ['rs', 'w m-2 h-1', None, 1, 0.0036],
    ]
)
def test_convert_solar_irradiance(variable, unit, timestep, value, expected):
    assert expected == units.convert(value, variable, unit, timestep)


@pytest.mark.parametrize(
    'variable, unit, value, expected',
    [
        ['uz', 'mph', 1, 0.44704],
    ]
)
def test_convert_velocity(variable, unit, value, expected):
    assert expected == units.convert(value, variable, unit)


@pytest.mark.parametrize(
    'variable, unit, value, expected',
    [
        ['zw', 'Feet', 1, 0.3048],
        ['elev', 'ft', 1, 0.3048],
    ]
)
def test_convert_length(variable, unit, value, expected):
    assert expected == units.convert(value, variable, unit)


@pytest.mark.parametrize(
    'variable, unit, value, expected',
    [
        ['lat', 'rad', math.pi, 180],
        ['lat', 'radian', math.pi, 180],
        ['lon', 'Radians', math.pi, 180],
    ]
)
def test_convert_angle(variable, unit, value, expected):
    assert expected == units.convert(value, variable, unit)

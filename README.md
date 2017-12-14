[![Build Status](https://travis-ci.org/Open-ET/RefET.svg?branch=master)](https://travis-ci.org/Open-ET/RefET)
[![PyPI version](https://badge.fury.io/py/RefET.svg)](https://badge.fury.io/py/RefET)

# ASCE Standardized Reference Evapotranspiration (ET)

NumPy functions for computing daily and hourly reference ET.

## Usage

#### Daily

To compute a single daily ETr value using weather data for 2015-07-01 from the [Fallon, NV AgriMet station](https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html).
The necessary unit conversions are shown on the input values.
The raw input data is available [here](https://www.usbr.gov/pn-bin/daily.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=ETRS&pcode=MN&pcode=MX&pcode=SR&pcode=YM&pcode=UA).

```
import math
import refet

# Unit conversions
tmin_c = (66.65 - 32) * (5.0 / 9)                          # F -> C
tmax_c = (102.80 - 32) * (5.0 / 9)                         # F -> C
tdew_c = (57.26 - 32) * (5.0 / 9)                          # F -> C
ea = 0.6108 * math.exp(17.27 * tdew_c / (tdew_c + 237.3))  # kPa
rs = (674.07 * 0.041868)                                   # Langleys -> MJ m-2 d-1
uz = 4.80 * 0.44704                                        # mpg -> m s-1
lat_radians = (39.4575 * math.pi / 180)                    # degrees -> radians

etr = refet.daily(
    tmin=tmin_c, tmax=tmax_c, ea=ea, rs=rs, uz=uz,
    zw=3, elev=1208.5, lat=lat_radians, doy=182, ref_type='etr')

print('ETr: {:.2f} mm'.format(float(etr)))
print('ETr: {:.2f} in'.format(float(etr) / 25.4))
```

#### Hourly

To compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the [Fallon, NV AgriMet station](https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html).
The necessary unit conversions are shown on the input values.
The raw input data is available [here](https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1).

```
import math
import refet

# Unit conversions
tmean_c = (91.80 - 32) * (5.0 / 9)           # F -> C
ea = 1.20                                    # kPa
rs = (61.16 * 0.041868)                      # W m-2 -> MJ m-2 h-1
uz = 3.33 * 0.44704                          # mph -> m s-1
lat_radians = (39.4575 * math.pi / 180)      # degrees -> radians
lon_radians = (-118.77388 * math.pi / 180)   # degrees -> radians

etr = refet.hourly(
    tmean=tmean_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
    lat=lat_radians, lon=lon_radians, doy=182, time=18, ref_type='etr')

print('ETr: {:.2f} mm'.format(float(etr)))
print('ETr: {:.2f} in'.format(float(etr) / 25.4))
```

## Limitations

The functions have **not** been tested for multi-dimensional arrays (i.e. time series or grids).

Currently the user must handle all of the file I/O and unit conversions.

### Units

The user must convert all input values to the expected units listed below:

type (variable) | units
-----|------
temperature (tmin, tmax, tmean) | C
actual vapor pressure (ea) | kPa
solar radiation (rs) | MJ m-2
wind speed (uz) | m s-1
wind speed measurement height (zw) | m
station elevation (elev) | m
latitude (lat) | radians
longitude (lon - hourly only) | radians
time (time - hourly only) | UTC hour at the start of the time period

### Cloudiness Fraction (hourly)

The hourly reference ET calculation is currently performed independently for each time step which causes the cloudiness fraction (fcd) calculation for very low sun angles to be incorrect.

## Installation

To install the RefET python module:
```
pip install refet
```

## Validation

Please see the [validation document](VALIDATION.md) for additional details on the source of the test values and the comparison of the functions to the Ref-ET software.

## Dependencies

* [numpy](http://www.numpy.org)

Modules needed to run the test suite:
* [pandas](http://pandas.pydata.org)
* [pytest](https://docs.pytest.org/en/latest/)
* [pytz](http://pythonhosted.org/pytz/)

## References

ASCE-EWRI Standardized Reference Evapotranspiration Equation (2005)
* [Report](http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf)
* [Appendix](http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf)

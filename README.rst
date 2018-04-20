===================================================
ASCE Standardized Reference Evapotranspiration (ET)
===================================================

|version| |build|

NumPy functions for computing daily and hourly reference ET.

Usage
=====

Daily Example
-------------

The following demonstrates how to compute a single daily ETr value using weather data for 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/daily.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=ETRS&pcode=MN&pcode=MX&pcode=SR&pcode=YM&pcode=UA>`__.

.. code-block:: console

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

    etr = refet.Daily(
        tmin=tmin_c, tmax=tmax_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat_radians, doy=182, method='asce').etr()

    print('ETr: {:.2f} mm'.format(float(etr)))

Hourly Example
--------------

The following demonstrates how to compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1>`__

.. code-block:: console

    import math
    import refet

    # Unit conversions
    tmean_c = (91.80 - 32) * (5.0 / 9)           # F -> C
    ea = 1.20                                    # kPa
    rs = (61.16 * 0.041868)                      # Langleys -> MJ m-2 h-1
    uz = 3.33 * 0.44704                          # mph -> m s-1
    lat_radians = (39.4575 * math.pi / 180)      # degrees -> radians
    lon_radians = (-118.77388 * math.pi / 180)   # degrees -> radians

    etr = refet.Hourly(
        tmean=tmean_c, ea=ea, rs=rs, uz=uz, zw=3, elev=1208.5,
        lat=lat_radians, lon=lon_radians, doy=182, time=18,
        method='asce').etr()

    print('ETr: {:.2f} mm'.format(float(etr)))


Input Parameters
================

Required Parameters (hourly & daily)
------------------------------------

========  ==========  ====================================================
Variable  Type        Description [units]
========  ==========  ====================================================
ea        ndarray     Actual vapor pressure [kPa]
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 day-1]
uz        ndarray     Wind speed [m/s]
zw        float       Wind speed height [m]
elev      ndarray     Elevation [m]
lat       ndarray     Latitude [radians]
doy       ndarray     Day of year
========  ==========  ====================================================

Required Daily Parameters
-------------------------

========  ==========  ====================================================
Variable  Type        Description [units]
========  ==========  ====================================================
tmin      ndarray     Minimum daily temperature [C]
tmax      ndarray     Maximum daily temperature [C]
========  ==========  ====================================================

Required Hourly Parameters
--------------------------

========  ==========  ====================================================
Variable  Type        Description [units]
========  ==========  ====================================================
tmean     ndarray     Average hourly temperature [C]
lon       ndarray     Longitude [radians]
time      ndarray     UTC hour at start of time period
========  ==========  ====================================================

Optional Parameters
-------------------

===========  ==========  ====================================================
Variable     Type        Description [units]
===========  ==========  ====================================================
method       str         | Calculation method

                         * 'asce' -- Calculations will follow ASCE-EWRI 2005 (default)
                         * 'refet' -- Calculations will follow RefET software

rso_type     str         | Override default clear sky solar radiation (Rso) calculation
                         | Defaults to None if not set

                         * 'full' -- Full clear sky solar formulation
                         * 'simple' -- Simplified clear sky solar formulation
                         * 'array' -- Read Rso values from "rso" function parameter

rso          array_like  | Clear sky solar radiation [MJ m-2 day-1]

                         * Only used if rso_type == 'array'
                         * Defaults to None if not set

input_units  dict        | Override default input unit types
                         | Input values will be converted to default unit types

===========  ==========  ====================================================

Issues
======

The functions have **not** been tested for inputs with different shapes/sizes and the broadcasting may not work correctly.

Currently the user must handle all of the file I/O and unit conversions.

The user must handle all QA/QC of the input data and no missing data will be filled.

Cloudiness Fraction (hourly)
----------------------------

The hourly reference ET calculation is currently performed independently for each time step which causes the cloudiness fraction (fcd) calculation for very low sun angles to be incorrect.

ASCE vs. RefET
==============

TODO Discuss differences between "asce" and "refet" methods.

Installation
============

To install the RefET python module:

.. code-block:: console

    pip install refet

Validation
==========

Please see the `validation document <VALIDATION.md>`__ for additional details on the source of the test values and the comparison of the functions to the Ref-ET software.

Dependencies
============

 * `numpy <http://www.numpy.org>`__

Modules needed to run the test suite:

 * `pandas <http://pandas.pydata.org>`__
 * `pytest <https://docs.pytest.org/en/latest/>`__
 * `pytz <http://pythonhosted.org/pytz/>`__

References
==========

ASCE-EWRI Standardized Reference Evapotranspiration Equation (2005)

 * `Report <http://www.kimberly.uidaho.edu/water/asceewri/ascestzdetmain2005.pdf>`__
 * `Appendix <http://www.kimberly.uidaho.edu/water/asceewri/appendix.pdf>`__

.. |build| image:: https://travis-ci.org/Open-ET/RefET.svg?branch=master
   :alt: Build status
   :target: https://travis-ci.org/Open-ET/RefET
.. |version| image:: https://badge.fury.io/py/RefET.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/RefET

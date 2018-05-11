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

    # Compute actual vapor pressure from Tdew
    tdew_c = (49.84 - 32) * (5.0 / 9)                          # F -> C
    ea = 0.6108 * math.exp(17.27 * tdew_c / (tdew_c + 237.3))  # kPa
    # ea = refet.calcs._saturated_vapor_pressure(tdew_c)

    etr = refet.Daily(
        tmin=66.65, tmax=102.80, ea=ea, rs=674.07, uz=4.80, zw=3, elev=1208.5,
        lat=39.4575, doy=182, method='asce',
        input_units={'tmin': 'F', 'tmax': 'F', 'rs': 'Langleys', 'uz': 'mph',
                     'lat': 'deg'}
        ).etr()

    print('ETr: {:.2f} mm'.format(float(etr)))

Hourly Example
--------------

The following demonstrates how to compute a single hourly ETr value using weather data for 18:00 UTC (11:00 AM PDT) on 2015-07-01 from the `Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are shown on the input values.
The raw input data is available `here <https://www.usbr.gov/pn-bin/instant.pl?station=FALN&year=2015&month=7&day=1&year=2015&month=7&day=1&pcode=OB&pcode=EA&pcode=WS&pcode=SI&print_hourly=1>`__

.. code-block:: console

    import refet

    etr = refet.Hourly(
        tmean=91.80, ea=1.20 , rs=61.16, uz=3.33, zw=3, elev=1208.5,
        lat=39.4575, lon=-118.77388, doy=182, time=18, method='asce',
        input_units={'tmean': 'F', 'rs': 'Langleys', 'uz': 'mph', 'lat': 'deg'}
        ).etr()

    print('ETr: {:.2f} mm'.format(float(etr)))


Input Parameters
================

Required Parameters (hourly & daily)
------------------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
ea        ndarray     Actual vapor pressure [kPa]
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 day-1]
uz        ndarray     Wind speed [m s-1]
zw        float       Wind speed height [m]
elev      ndarray     Elevation [m]
lat       ndarray     Latitude [degrees]
doy       ndarray     Day of year
========  ==========  ====================================================

Required Daily Parameters
-------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
tmin      ndarray     Minimum daily temperature [C]
tmax      ndarray     Maximum daily temperature [C]
========  ==========  ====================================================

Required Hourly Parameters
--------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
tmean     ndarray     Average hourly temperature [C]
lon       ndarray     Longitude [degrees]
time      ndarray     UTC hour at start of time period
========  ==========  ====================================================

Optional Parameters
-------------------

===========  ==========  ====================================================
Variable     Type        Description [default units]
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

The cloudiness fraction (fcd) is computed as the ratio of the measured solar radiation (Rs) to the theoretical clear sky solar radiation (Rso).  This ratio cannot be computed directly at night since Rso is 0.  ASCE-EWRI 2005 suggests computing a representative nighttime fcd based on the fcd at sunset and/or sunrise.

In the RefET module fcd is hard coded to 1 for all time steps with very low sun angles since the hourly reference ET is computed independently for each time step.

ASCE vs. RefET
==============

TODO: Discuss the differences between "asce" and "refet" methods.

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

.. |build| image:: https://travis-ci.org/DRI-WSWUP/RefET.svg?branch=master
   :alt: Build status
   :target: https://travis-ci.org/DRI-WSWUP/RefET
.. |version| image:: https://badge.fury.io/py/RefET.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/RefET

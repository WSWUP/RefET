===================================================
ASCE Standardized Reference Evapotranspiration (ET)
===================================================

|version| |build|

NumPy functions for computing daily and hourly reference ET following the ASCE Standardized Reference Evapotranspiration Equations (ASCE2005_).

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

    # The actual vapor pressure could be computed from the dew point temperature below
    #   or the tdew can be passed directly to the function
    # Convert the dew point temperature to Celsius
    # tdew = units._f2c(49.84)
    # ea = 0.6108 * math.exp(17.27 * tdew / (tdew + 237.3))
    # ea = refet.calcs._sat_vapor_pressure(tdew)

    etr = refet.Daily(
        tmin=66.65, tmax=102.80, tdew=49.84, rs=674.07, uz=4.80,
        zw=3, elev=1208.5, lat=39.4575, doy=182, method='asce',
        input_units={'tmin': 'F', 'tmax': 'F', 'tdew': 'F', 'rs': 'Langleys',
                     'uz': 'mph', 'lat': 'deg'}
        ).etr()

    print(f'ETr: {float(etr):.2f} mm')

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

    print(f'ETr: {float(etr):.2f} mm')


Input Parameters
================

Required Parameters (hourly & daily)
------------------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
uz        ndarray     Wind speed [m s-1]
zw        float       Wind speed height [m]
elev      ndarray     Elevation [m]
lat       ndarray     Latitude [degrees]
doy       ndarray     Day of year
========  ==========  ====================================================

Required Ea Parameters (hourly & daily)
---------------------------------------------------

Either the "ea" or "tdew" parameter must be set

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
ea        ndarray     Actual vapor pressure [kPa]
tdew      ndarray     Dew point temperature [C]
========  ==========  ====================================================

Required Daily Parameters
-------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 d-1]
tmin      ndarray     Minimum daily temperature [C]
tmax      ndarray     Maximum daily temperature [C]
========  ==========  ====================================================

Required Hourly Parameters
--------------------------

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 h-1]
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

rso          array_like  | Clear sky solar radiation [MJ m-2 d-1 or MJ m-2 h-1]

                         * Only used if rso_type == 'array'
                         * Defaults to None if not set

input_units  dict        | Override default input unit types
                         | Input values will be converted to default unit types

===========  ==========  ====================================================

Installation
============

The RefET python module can be installed with conda or pip:

.. code-block:: console

    conda install refet

Issues
======

The functions have **not** been tested for inputs with different shapes/sizes and the broadcasting may not work correctly.

The user must handle the following:
 + File I/O
 + QA/QC of the input data
 + Filling missing or bad data

Cloudiness Fraction (hourly)
----------------------------

The cloudiness fraction (fcd) is computed as the ratio of the measured solar radiation (Rs) to the theoretical clear sky solar radiation (Rso).  This ratio cannot be computed directly at night since Rso is 0.  ASCE2005_ suggests computing a representative nighttime fcd based on the fcd at sunset and/or sunrise.

In the RefET module fcd is hard coded to 1 for all time steps with very low sun angles since the hourly reference ET is computed independently for each time step.

Calculation Method - ASCE vs. RefET
===================================

The main difference between the two "methods" is that the "asce" method attempts to follow the equations in ASCE2005_, whereas the "refet" method attempts to follow the calculations of the `RefET Software <https://www.uidaho.edu/cals/kimberly-research-and-extension-center/research/water-resources/ref-et-software>`__ as closely as possible.  The difference in output between these methods is generally negligible (if not identical for realistic numbers of significant digits).  Note that the default is set to "asce" to best match the calculations a user would expect to have happen. The "refet" method was added in order to help validate this code to the RefET Software.

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

.. _references:

.. [ASCE2005]
 | ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration equation.
 | `https://ascelibrary.org/doi/book/10.1061/9780784408056 <https://ascelibrary.org/doi/book/10.1061/9780784408056>`__

.. |build| image:: https://github.com/WSWUP/refet/workflows/build/badge.svg
   :alt: Build status
   :target: https://github.com/WSWUP/refet
.. |version| image:: https://badge.fury.io/py/refet.svg
   :alt: Latest version on PyPI
   :target: https://badge.fury.io/py/refet

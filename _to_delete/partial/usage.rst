Usage
=====

Daily example
-------------

Compute a single daily ETr value using weather data for 2015-07-01 from the
`Fallon, NV AgriMet station <https://www.usbr.gov/pn/agrimet/agrimetmap/falnda.html>`__.
The necessary unit conversions are handled through ``input_units``.

.. code-block:: python

    import refet

    etr = refet.Daily(
        tmin=66.65, tmax=102.80, tdew=49.84, rs=674.07, uz=4.80,
        zw=3, elev=1208.5, lat=39.4575, doy=182, method='asce',
        input_units={'tmin': 'F', 'tmax': 'F', 'tdew': 'F', 'rs': 'Langleys',
                     'uz': 'mph', 'lat': 'deg'}
        ).etr()

    print(f'ETr: {float(etr):.2f} mm')

Hourly example
--------------

Compute a single hourly ETr value for 18:00 UTC (11:00 AM PDT) on 2015-07-01
at the same station.

.. code-block:: python

    import refet

    etr = refet.Hourly(
        tmean=91.80, ea=1.20, rs=61.16, uz=3.33, zw=3, elev=1208.5,
        lat=39.4575, lon=-118.77388, doy=182, time=18, method='asce',
        input_units={'tmean': 'F', 'rs': 'Langleys', 'uz': 'mph', 'lat': 'deg'}
        ).etr()

    print(f'ETr: {float(etr):.2f} mm')

Intermediate terms
------------------

The ``results()`` method returns the input and intermediate terms computed by
the constructor (net radiation, vapor pressure deficit, clear sky solar
radiation, and so on) as a dict of arrays, which is useful for diagnostics and
plotting:

.. code-block:: python

    d = refet.Daily(
        tmin=19.25, tmax=39.33, tdew=9.91, rs=28.22, uz=2.15,
        zw=3, elev=1208.5, lat=39.4575, doy=182)
    terms = d.results()
    print(sorted(terms.keys()))

Input parameters
----------------

Required parameters (hourly and daily)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
uz        ndarray     Wind speed [m s-1]
zw        float       Wind speed height [m]
elev      ndarray     Elevation [m]
lat       ndarray     Latitude [degrees]
doy       ndarray     Day of year
========  ==========  ====================================================

Required vapor pressure parameters (hourly and daily)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Either the "ea" or "tdew" parameter must be set.

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
ea        ndarray     Actual vapor pressure [kPa]
tdew      ndarray     Dew point temperature [C]
========  ==========  ====================================================

Required daily parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 d-1]
tmin      ndarray     Minimum daily temperature [C]
tmax      ndarray     Maximum daily temperature [C]
========  ==========  ====================================================

Required hourly parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

========  ==========  ====================================================
Variable  Type        Description [default units]
========  ==========  ====================================================
rs        ndarray     Incoming shortwave solar radiation [MJ m-2 h-1]
tmean     ndarray     Average hourly temperature [C]
lon       ndarray     Longitude [degrees]
time      ndarray     UTC hour at start of time period
========  ==========  ====================================================

Optional parameters
~~~~~~~~~~~~~~~~~~~

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

Limitations
-----------

The user must handle file I/O, quality control of the input data, and filling
of missing or bad values. Broadcasting of mixed input shapes is not fully
tested.

The hourly cloudiness fraction (fcd) is hard coded to 1 for time steps with
very low sun angles, since the hourly reference ET is computed independently
for each time step (ASCE-EWRI 2005 instead suggests carrying a representative
nighttime fcd from conditions at sunset or sunrise).

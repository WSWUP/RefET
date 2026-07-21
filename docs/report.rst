Station reports and AgriMet data
================================

The ``refet[report]`` extra adds two things: a client for the Bureau of
Reclamation AgriMet network (``refet.io.agrimet``) and a report builder
(``refet.report``) that turns a small TOML config into a self-contained,
interactive HTML report of daily and hourly reference ET.

.. code-block:: console

    pip install refet[report]
    refet-report stations/faln_2024_live.toml

A minimal config for a live AgriMet station:

.. code-block:: toml

    start = 2024-04-01
    end = 2024-09-30

    [source]
    type = "agrimet"
    station = "FALN"
    cache = ".cache"

    [daily.columns]
    tmin = "MN"
    tmax = "MX"
    tdew = "YM"
    uz = "UA"
    rs = "SR"

    [hourly.columns]
    tmean = "OB"
    tdew = "TP"
    uz = "WS"
    rs = "SI"

    [units]
    tmin = "F"
    tmax = "F"
    tmean = "F"
    tdew = "F"
    uz = "mph"
    rs = "Langleys"

Latitude, longitude, elevation, timezone, and the display name come from the
published AgriMet station list. Anemometer height defaults to the network
standard of 3 m; override it in ``[station] zw`` if a site differs. Local CSV
sources and validation against Ref-ET ``.out`` files are also supported; see
the worked examples in the repository's ``examples/stations/`` directory.

Every number in the report's prose is a derived fact computed from the data,
and the build fails rather than rendering a page whose narrative and numbers
could disagree. Reports for stations without reference data automatically
omit the validation panel.

The AgriMet client can also be used directly:

.. code-block:: python

    import datetime as dt
    from refet.io import agrimet

    daily = agrimet.fetch_daily(
        'FALN', dt.date(2015, 7, 1), dt.date(2015, 7, 31),
        ['MN', 'MX', 'SR', 'YM', 'UA'], cache='.cache')

Values come back in native AgriMet units (F, Langleys, mph) with every known
missing-data convention normalized to NaN, ready to feed
:class:`refet.Daily` with the matching ``input_units`` mapping.

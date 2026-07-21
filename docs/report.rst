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

Config reference
----------------

Relative paths in a config resolve against the config file's directory.
Top-level keys:

=============  ==========================================================
Key            Meaning
=============  ==========================================================
start, end     Inclusive date range as TOML dates (required)
title          Page title (auto-filled for AgriMet sources)
name           Station display name (auto-filled for AgriMet sources)
network        Network label shown in the report (auto-filled)
timezone       IANA zone of the raw timestamps, the DST-observing zone,
               not a fixed offset (auto-filled for AgriMet sources)
method         'asce' (default) or 'refet', passed to the compute classes
rso_type       Clear sky solar model override: 'simple' or 'full'
               (defaults to the method's standard choice)
=============  ==========================================================

Sections:

==================  =====================================================
Section             Meaning
==================  =====================================================
[source]            type = "csv" (default) or "agrimet". AgriMet sources
                    need a station id and accept an optional cache
                    directory.
[station]           lat, lon, elev [m], zw [m] overrides. Required for
                    csv sources with no reference file to scrape them
                    from; optional otherwise (explicit values always win).
[daily]             csv path (csv sources) and an optional Ref-ET .out
                    reference file; with a reference, the report gains a
                    validation panel.
[daily.columns]     report variable = source column name. Must map tmin,
                    tmax, uz, rs, and either tdew or ea.
[hourly]            As [daily], for the hourly source.
[hourly.columns]    Must map tmean, uz, rs, and either tdew or ea.
[units]             Unit string per variable, anything
                    refet.units.convert accepts. Omit for library
                    defaults (C, m s-1, MJ m-2).
==================  =====================================================

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

RefET
=====

NumPy functions for computing daily and hourly reference evapotranspiration
(ET) following the ASCE Standardized Reference Evapotranspiration Equations
(ASCE-EWRI 2005).

The package computes standardized reference ET for the short (grass, ETo) and
tall (alfalfa, ETr) reference surfaces from standard weather station inputs,
and is validated against the University of Idaho Ref-ET software using a full
year of AgriMet data (see :doc:`validation`).

Beyond the core equations, optional modules cover the rest of the station
workflow: an AgriMet data client and a self-contained HTML report builder
(:doc:`report`), input screening flags (``refet.qaqc``), and documented
estimators for missing inputs (``refet.estimate``).

For the equivalent functions implemented for Google Earth Engine, see
`openet-refet-gee <https://github.com/Open-ET/openet-refet-gee>`__.

.. toctree::
   :maxdepth: 1

   installation
   usage
   report
   api
   fao56
   validation
   development

References
----------

ASCE-EWRI (2005). The ASCE standardized reference evapotranspiration equation.
https://ascelibrary.org/doi/book/10.1061/9780784408056

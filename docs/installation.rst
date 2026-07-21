Installation
============

RefET can be installed with pip:

.. code-block:: console

    pip install refet

or from conda-forge:

.. code-block:: console

    conda install conda-forge::refet

The core package depends only on `NumPy <https://numpy.org>`__ (>=2.0) and
requires Python 3.11 or newer.

Optional extras
---------------

The data access and reporting modules have extra dependencies, installed on
demand so the core stays lightweight:

``refet[data]``
    The ``refet.io.agrimet`` AgriMet client (adds pandas).

``refet[report]``
    The ``refet.report`` builder and ``refet-report`` command (adds pandas
    and tzdata).

``refet[test]``
    Everything the test suite needs (pandas, pytz, pytest, tzdata).

``refet[docs]``
    The documentation build (sphinx, furo, myst-parser, pandas).

The pure-NumPy modules (``refet.calcs``, ``refet.units``, ``refet.qaqc``,
``refet.estimate``) work with the base install.

For a development install, see :doc:`development`.

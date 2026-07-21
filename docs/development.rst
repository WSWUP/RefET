Development
===========

Working on the package
----------------------

.. code-block:: console

    git clone https://github.com/WSWUP/RefET.git
    cd RefET
    pip install -e .[test]
    python -m pytest

The suite runs entirely offline against the recorded Fallon 2015 data. One
test fetches live AgriMet data and is skipped unless you opt in:

.. code-block:: console

    RUN_NETWORK_TESTS=1 python -m pytest tests/test_io_agrimet.py

Lint matches CI:

.. code-block:: console

    pip install ruff
    ruff check .

Docs build locally with:

.. code-block:: console

    pip install -e .[docs]
    sphinx-build -b html docs docs/_build/html

The science guardrail
---------------------

``refet.calcs``, ``refet.daily``, ``refet.hourly``, and ``refet.units`` are
the validated core. The full-year Fallon 2015 golden tests
(``tests/test_daily_refet_output.py`` and ``tests/test_hourly_refet_output.py``)
compare every computed value against University of Idaho Ref-ET software
output; any change that shifts those results needs an explicit, documented
reason (the CHANGELOG entry for the hourly VPD clamp in 0.6.0 is the model).
New capability belongs in separate modules with optional dependencies, not
in the core.

Continuous integration
----------------------

Three workflows run from ``.github/workflows/``:

- **tests.yml**: ruff lint plus the pytest suite on Python 3.11 through 3.14
  across Ubuntu and Windows, on every push and pull request to main.
- **docs.yml**: builds this documentation and deploys it to GitHub Pages on
  pushes to main. One-time repository setup: Settings, Pages, set Source to
  "GitHub Actions".
- **publish.yml**: builds the sdist and wheel and uploads to PyPI when a
  GitHub release is published, using PyPI Trusted Publishing. One-time setup
  on PyPI: project Publishing settings, add a trusted publisher for this
  repository with workflow ``publish.yml`` and environment ``pypi``.

Cutting a release
-----------------

1. Update the version in ``pyproject.toml`` and ``CITATION.cff``, date the
   release section in ``CHANGELOG.md``, and merge to main.
2. Confirm CI is green and the docs deployed.
3. Create a GitHub release with a ``v`` tag (for example ``v0.6.0``); the
   publish workflow uploads to PyPI automatically.
4. The conda-forge feedstock picks up the new PyPI release through its bot;
   review the automated PR it opens.

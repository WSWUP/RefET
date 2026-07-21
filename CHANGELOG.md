# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [semantic versioning](https://semver.org/).

## [0.6.0] - Unreleased

### Fixed
- `Daily` with `rso_type='simple'` now uses the unit-converted elevation.
  Previously the simple Rso branch read the raw `elev` argument, so an
  elevation supplied in feet via `input_units` skipped conversion.
- `Daily` with `rso_type='array'` converts the supplied `rso` to an ndarray
  and raises `ValueError` when `rso` is missing instead of computing with
  `None`.
- `Hourly` no longer discards the ndarray conversion of `doy`.
- `units.convert()` no longer modifies caller arrays in place, and integer
  inputs no longer raise a casting error.
- Two shadowed duplicate test functions (`test_refet_daily_tdew`,
  `test_refet_hourly_tdew`) were renamed so both tests actually run.

### Changed
- `Hourly` vapor pressure deficit is clamped at zero, matching `Daily` and
  ASCE-EWRI (2005) intent. This only affects time steps where Tdew exceeds
  Tmean (a sensor artifact). In the Fallon 2015 validation year this touches
  164 of 8,758 hours, all at night, none in the Ref-ET golden comparison set.
- `Daily` and `Hourly` raise `ValueError` instead of a bare `Exception` when
  neither `ea` nor `tdew` is given.
- The `input_units` default changed from a mutable `{}` to `None`.
- Tests import the installed `refet` package instead of the source tree.
- CI runs on Python 3.11 through 3.14 on Ubuntu and Windows, with ruff lint.
- Packaging metadata modernized: SPDX license expression (PEP 639),
  version classifiers, keywords.

### Added
- `results()` method on `Daily` and `Hourly` returning the input and
  intermediate terms (net radiation, VPD, clear sky solar, and so on) as a
  dict of arrays.
- Type hints on the public API and a `py.typed` marker (PEP 561), so editors
  and type checkers see real signatures.
- Sphinx documentation site (usage, API reference, validation, and a note on
  the relationship between ASCE standardized reference ET and FAO-56) with a
  GitHub Pages deploy workflow.
- `CITATION.cff` citation metadata and this changelog.
- PyPI publish workflow using trusted publishing, triggered by GitHub
  releases.
- `__all__` in the package init.

## [0.5.0] - 2026-05-20 (not published to PyPI)

### Changed
- Minimum NumPy version raised to 2.0.
- Minimum Python version raised to 3.11.
- `__version__` is read from the installed package metadata.
- Requirements files removed in favor of `pyproject.toml`.

## [0.4.2] - 2023-06-30

Last release published to PyPI and conda-forge before this changelog existed.

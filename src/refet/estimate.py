"""Estimated inputs for gap filling, clearly separated from measurements.

Everything in this module is an estimate with documented assumptions, per
FAO-56 and ASCE-EWRI (2005) Appendix E. Use these to fill gaps or run
stations with missing sensors, and say so in your methods section; never
mix them silently with measured values.
"""
import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def tdew_from_tmin(tmin: ArrayLike, ko: float = 2.0) -> FloatArray | float:
    """Estimate dew point from minimum air temperature.

    Parameters
    ----------
    tmin : scalar or array_like of shape(M, )
        Minimum daily air temperature [C].
    ko : float, optional
        Offset [C] subtracted from tmin. ASCE-EWRI (2005) Appendix E
        suggests 2 to 4 C for arid and semiarid sites and 0 for humid
        sites; FAO-56 (Annex 6) uses the same approach.

    Returns
    -------
    ndarray
        Estimated dew point temperature [C].

    """
    return np.asarray(tmin, dtype=float) - ko


def tdew_from_ea(ea: ArrayLike) -> FloatArray | float:
    """Dew point from actual vapor pressure, inverting the Magnus form.

    Parameters
    ----------
    ea : scalar or array_like of shape(M, )
        Actual vapor pressure [kPa].

    Returns
    -------
    ndarray
        Dew point temperature [C].

    Notes
    -----
    Inverts es = 0.6108 * exp(17.27 * T / (T + 237.3)) (FAO-56 Eq. 3-11).
    This is a unit conversion more than an estimate, but it lives here so
    the inverse relationship stays next to its assumptions.

    """
    x = np.log(np.asarray(ea, dtype=float) / 0.6108)
    return 237.3 * x / (17.27 - x)


def rs_hargreaves(ra: ArrayLike, tmin: ArrayLike, tmax: ArrayLike,
                  krs: float = 0.16) -> FloatArray | float:
    """Estimate solar radiation from the daily temperature range.

    Hargreaves-Samani radiation formula (FAO-56 Eq. 50; ASCE-EWRI 2005
    Appendix E): rs = krs * sqrt(tmax - tmin) * ra.

    Parameters
    ----------
    ra : scalar or array_like of shape(M, )
        Extraterrestrial radiation [MJ m-2 d-1], e.g. from
        :func:`refet.calcs.ra_daily`.
    tmin : scalar or array_like of shape(M, )
        Minimum daily air temperature [C].
    tmax : scalar or array_like of shape(M, )
        Maximum daily air temperature [C].
    krs : float, optional
        Adjustment coefficient: 0.16 for interior locations, 0.19 for
        coastal locations.

    Returns
    -------
    ndarray
        Estimated incoming solar radiation [MJ m-2 d-1].

    """
    dt = np.asarray(tmax, dtype=float) - np.asarray(tmin, dtype=float)
    return krs * np.sqrt(np.maximum(dt, 0)) * np.asarray(ra, dtype=float)

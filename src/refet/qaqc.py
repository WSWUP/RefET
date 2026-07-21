"""Input data screening flags, following ASCE-EWRI (2005) Appendix D.

These functions flag suspect inputs; they never modify data. Filling,
correction, and drift analysis are deliberately out of scope (see WSWUP's
agweather-qaqc for a full correction workflow). The typical use is a check
before computing ET from station data you have not looked at::

    d = refet.Daily(...)
    flags = refet.qaqc.check(d)
    for name, flag in flags.items():
        if flag.any():
            print(f'{name}: {int(flag.sum())} flagged')

Every function returns a boolean array shaped like its inputs, True where
the value is suspect.
"""
import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import calcs

BoolArray = NDArray[np.bool_]

# Measured Rs above this multiple of clear sky Rso indicates a calibration
# or timing problem (the envelope check of ASCE-EWRI 2005 Appendix D).
RS_RSO_TOL = 1.1
# Wind speeds above this [m s-1] are treated as sensor spikes.
UZ_MAX = 30.0


def flag_rs_above_rso(rs: ArrayLike, rso: ArrayLike,
                      tol: float = RS_RSO_TOL) -> BoolArray:
    """Measured solar above tol * clear sky solar (envelope check)."""
    return np.asarray(rs, dtype=float) > tol * np.asarray(rso, dtype=float)


def flag_negative(values: ArrayLike) -> BoolArray:
    """Values below zero, for variables that cannot be negative."""
    return np.asarray(values, dtype=float) < 0


def flag_range(values: ArrayLike, lo: float, hi: float) -> BoolArray:
    """Values outside [lo, hi]."""
    v = np.asarray(values, dtype=float)
    return (v < lo) | (v > hi)


def flag_tmax_not_above_tmin(tmax: ArrayLike, tmin: ArrayLike) -> BoolArray:
    """Daily maximum temperature at or below the minimum."""
    return np.asarray(tmax, dtype=float) <= np.asarray(tmin, dtype=float)


def flag_ea_above_saturation(ea: ArrayLike, temperature: ArrayLike) -> BoolArray:
    """Vapor pressure above saturation at the given temperature.

    For daily data pass tmax (ea should stay below es(tmax)); for hourly
    pass tmean. True flags indicate a humidity or temperature sensor
    problem; the ET classes clamp the resulting negative VPD to zero, so
    this is where those time steps become visible.
    """
    return np.asarray(ea, dtype=float) > calcs.sat_vapor_pressure(temperature)


def check(model, rs_rso_tol: float = RS_RSO_TOL,
          uz_max: float = UZ_MAX) -> dict[str, BoolArray]:
    """Screen the inputs of a Daily or Hourly instance.

    Parameters
    ----------
    model : refet.Daily or refet.Hourly
        A constructed instance; its converted inputs and computed Rso are
        read directly.
    rs_rso_tol : float, optional
        Tolerance for the Rs vs Rso envelope check.
    uz_max : float, optional
        Maximum plausible wind speed [m s-1].

    Returns
    -------
    dict of ndarray
        Boolean flag arrays keyed by check name. True marks a suspect
        time step. No data is modified.

    """
    flags = {}
    daily = hasattr(model, 'tmax')
    if daily:
        flags['tmax_not_above_tmin'] = flag_tmax_not_above_tmin(
            model.tmax, model.tmin)
        flags['ea_above_saturation'] = flag_ea_above_saturation(
            model.ea, model.tmax)
    else:
        flags['ea_above_saturation'] = flag_ea_above_saturation(
            model.ea, model.tmean)
    flags['rs_negative'] = flag_negative(model.rs)
    flags['rs_above_rso'] = flag_rs_above_rso(model.rs, model.rso, rs_rso_tol)
    flags['ea_negative'] = flag_negative(model.ea)
    flags['uz_negative'] = flag_negative(model.uz)
    flags['uz_above_max'] = np.asarray(model.uz, dtype=float) > uz_max
    return flags


def counts(flags: dict[str, BoolArray]) -> dict[str, int]:
    """Number of flagged time steps per check, for a quick summary."""
    return {name: int(np.asarray(flag).sum()) for name, flag in flags.items()}

import numpy as np
import pytest

from refet import calcs, estimate


def test_tdew_from_tmin():
    assert estimate.tdew_from_tmin(10.0) == 8.0
    assert estimate.tdew_from_tmin(10.0, ko=0.0) == 10.0
    out = estimate.tdew_from_tmin(np.array([5.0, 15.0]), ko=3.0)
    assert out.tolist() == [2.0, 12.0]


@pytest.mark.parametrize('t', [-20.0, 0.0, 10.0, 25.0, 40.0])
def test_tdew_from_ea_round_trip(t):
    """Inverting the saturation curve recovers the temperature"""
    ea = float(calcs.sat_vapor_pressure(t)[0])
    assert estimate.tdew_from_ea(ea) == pytest.approx(t, abs=1e-6)


def test_rs_hargreaves():
    # FAO-56 Eq. 50 with a 20 C range: 0.16 * sqrt(20) * ra
    ra = 40.0
    expected = 0.16 * np.sqrt(20.0) * ra
    assert float(estimate.rs_hargreaves(ra, 10.0, 30.0)) == pytest.approx(expected)
    # Coastal coefficient
    assert float(estimate.rs_hargreaves(ra, 10.0, 30.0, krs=0.19)) == \
        pytest.approx(0.19 * np.sqrt(20.0) * ra)
    # Inverted temperatures do not produce NaN
    assert float(estimate.rs_hargreaves(ra, 30.0, 10.0)) == 0.0


def test_rs_hargreaves_vs_measured_scale():
    """Sanity: the estimate lands in the right range for a clear July day"""
    ra = float(np.asarray(calcs.ra_daily(39.4575 * np.pi / 180, 182)).ravel()[0])
    rs_est = float(estimate.rs_hargreaves(ra, 19.25, 39.33))
    # Measured Rs at Fallon on 2015-07-01 was ~28.2 MJ m-2 d-1
    assert 20.0 < rs_est < 35.0

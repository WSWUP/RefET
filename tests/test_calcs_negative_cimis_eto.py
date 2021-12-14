import math

import pytest

import refet.calcs as calcs

# CIMIS Data for 2021-12-05 at -120.12235, 36.33455
# Tdew: 8.329306602478027
# Tx: 9.572299003601074
# Tn: 7.329018592834473
# Rnl: 5.315489768981934
# Rs: -2.1786582469940186
# K: -0.40143924951553345
# U2: 1.3690364360809326
# ETo: 0.664989709854126
# ETo_ASCE: -0.33529889583587646
# ETr_ASCE: -0.3130447566509247
# Rso: (-2.1786582469940186 / -0.40143924951553345) = 5.427118174476651


@pytest.mark.parametrize(
    'lat, doy, ra',
    [
        # This is the expected Ra for the lat and doy
        [36.333404541015625 * (math.pi / 180), 339, 16.29900171668147],
        # Check the value during the summer
        [36.333404541015625 * (math.pi / 180), 200, 40.503108049470946],
    ])
def test_ra_daily(lat, doy, ra):
    assert float(calcs._ra_daily(lat, doy, method='asce')) == pytest.approx(ra)


@pytest.mark.parametrize(
    'ra, elev, rso',
    [
        # This is the expected Rso for the elev (and lat, doy)
        [16.29900171668147, 86.375, 12.252407812976669],
        # # This is the Ra needed to get the Rso computed from K and Rs
        # [7.22, 86.375, 5.427118174476651],
        # Check the value during the summer
        [40.503108049470946, 86.375, 30.44730015625867],
    ])
def test_rso_simple(ra, elev, rso):
    assert float(calcs._rso_simple(ra, elev)) == pytest.approx(rso)


@pytest.mark.parametrize(
    'rs, rso, fcd',
    [
        # Use the expected Rso based on the lat, doy, elev
        [-2.1786582469940186, 12.252407812976669, 0.055],
        # Use the Rso in the CIMIS inputs files
        [-2.1786582469940186, 5.427118174476651, 0.055],

        # Check that an Rs of 0 is clamped
        [0, 5.5, 0.055],
        # Check the summer
        [31, 30.44730015625867, 1.0],
    ])
def test_fcd_daily(rs, rso, fcd):
    assert float(calcs._fcd_daily(rs, rso)) == pytest.approx(fcd)


@pytest.mark.parametrize(
    'tmin, tmax, ea, fcd, rnl, method',
    [
        # This is the Rnl that is generated using the CIMIS inputs
        [7.329018592834473, 9.572299003601074,
         calcs._sat_vapor_pressure(8.329306602478027), 0.055,
         0.3278366577950547, 'asce'],
        # The fcd is not clamped in Spatial CIMIS calculation,
        #   so recompute it so that it can be negative
        # The fcd or Rnl sign needs to be flipped to match the values in the CIMIS files
        #   i.e. Rnl is positive in the file
        [7.329018592834473, 9.572299003601074,
         calcs._sat_vapor_pressure(8.329306602478027),
         0.8919429868459702,
         # -1 * (1.35 * -0.40143924951553345 - 0.35),
         # -1 * (1.35 * (-2.1786582469940186 / 5.427118174476651) - 0.35),
         5.315489768981934, 'cimis'],
        # The "asce" rnl calculation would be just a little bit different
        [7.329018592834473, 9.572299003601074,
         calcs._sat_vapor_pressure(8.329306602478027),
         0.8919429868459702, 5.316574686387661, 'asce'],

        # # This is the fcd value needed to get the Rnl on the CIMIS site
        # [7.329018592834473, 9.572299003601074,
        #  calcs._sat_vapor_pressure(8.329306602478027), 0.872, 5.315489768981934],
        # # Try using a negatice fcd
        # [7.329018592834473, 9.572299003601074,
        #  calcs._sat_vapor_pressure(8.329306602478027), -0.5, -3.046932843518462],
        # Testing out what the Rnl would be for an fcd of 1
        # [7.329018592834473, 9.572299003601074,
        #  calcs._sat_vapor_pressure(8.329306602478027), 1, 6.093865687036924],

     ])
def test_rnl_daily(tmin, tmax, ea, fcd, rnl, method):
    assert float(calcs._rnl_daily(tmax, tmin, ea, fcd, method)) == pytest.approx(rnl)


@pytest.mark.parametrize(
    'rs, rnl, rn',
    [
        # Rn generated using refet and the CIMIS inputs
        [-2.1786582469940186, 0.3278366577950547, -2.005403507980449],
        # Rn value used in the spatial CIMIS calculation for ETo
        #   (sign convention for Rnl in equation is flipped)
        [-2.1786582469940186, -1 * 5.315489768981934, 3.637922918796539],
        # Try clamping Rs to >= 0 to see if RefET value would be valid
        [0, 0.3278366577950547, -0.3278366577950547],
    ])
def test_rn_daily_negative_rs(rs, rnl, rn):
    assert float(calcs._rn_daily(rs, rnl)) == pytest.approx(rn)


 # This is the RefET computed value: -0.33529889583587646
 # This is the Spatial CIMIS value: 0.664989709854126
@pytest.mark.parametrize(
    'rn, g, tmean, u2, es, ea, es_slope, pair, cn, cd, etsz',
    [
        # This is the Rn using the Rs and Rnl from CIMIS input values
        [calcs._rn_daily(rs=-2.1786582469940186,
                         rnl=calcs._rnl_daily(
                             tmax=9.572299003601074, tmin=7.329018592834473,
                             ea=calcs._sat_vapor_pressure(8.329306602478027),
                             fcd=0.055, method='asce')),
         0, 0.5 * (7.329018592834473 + 9.572299003601074), 1.3690364360809326,
         0.5 * (calcs._sat_vapor_pressure(7.329018592834473) +
                calcs._sat_vapor_pressure(9.572299003601074)),
         calcs._sat_vapor_pressure(8.329306602478027),
         calcs._es_slope(0.5 * (7.329018592834473 + 9.572299003601074), method='asce'),
         calcs._air_pressure(86.375, method='asce'), 900, 0.34,
         -0.3353167433284481],
         # -0.33529889583587646],

        # This is the Rn to get the precomputed CIMIS ETo value
        [calcs._rn_daily(rs=-2.1786582469940186,
                         rnl=-1 * calcs._rnl_daily(
                             tmax=9.572299003601074, tmin=7.329018592834473,
                             ea=calcs._sat_vapor_pressure(8.329306602478027),
                             fcd=0.8919429868459702, method='cimis')),
         0, 0.5 * (7.329018592834473 + 9.572299003601074), 1.3690364360809326,
         0.5 * (calcs._sat_vapor_pressure(7.329018592834473) +
                calcs._sat_vapor_pressure(9.572299003601074)),
         calcs._sat_vapor_pressure(8.329306602478027),
         calcs._es_slope(0.5 * (7.329018592834473 + 9.572299003601074), method='cimis'),
         calcs._air_pressure(86.375, method='asce'), 900, 0.34,
         0.664871664059056],
         # 0.664989709854126],

        # Clamping the Rn >= 0 returns a positive ET
        [0, 0, 0.5 * (7.329018592834473 + 9.572299003601074), 1.3690364360809326,
         0.5 * (calcs._sat_vapor_pressure(7.329018592834473) +
                calcs._sat_vapor_pressure(9.572299003601074)),
         calcs._sat_vapor_pressure(8.329306602478027),
         calcs._es_slope(0.5 * (7.329018592834473 + 9.572299003601074), method='asce'),
         calcs._air_pressure(86.375, method='asce'), 900, 0.34, 0.0201],

        # Rn is still negative even if Rs is clamped to >= 0 (Rs computed above)
        [calcs._rn_daily(rs=0,
                         rnl=calcs._rnl_daily(
                             tmax=9.572299003601074, tmin=7.329018592834473,
                             ea=calcs._sat_vapor_pressure(8.329306602478027),
                             fcd=0.055, method='asce')),
         0, 0.5 * (7.329018592834473 + 9.572299003601074), 1.3690364360809326,
         0.5 * (calcs._sat_vapor_pressure(7.329018592834473) +
                calcs._sat_vapor_pressure(9.572299003601074)),
         calcs._sat_vapor_pressure(8.329306602478027),
         calcs._es_slope(0.5 * (7.329018592834473 + 9.572299003601074), method='asce'),
         calcs._air_pressure(86.375, method='asce'), 900, 0.34, -0.03800199489054891],
     ])
def test_etsz(rn, g, tmean, u2, es, ea, es_slope, pair, cn, cd, etsz):
    output = calcs._etsz(rn, g, tmean, u2, es - ea, es_slope, 0.000665 * pair, cn, cd)
    assert float(output) == pytest.approx(etsz, abs=0.0001)

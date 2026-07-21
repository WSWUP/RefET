import numpy as np
import pytest

import refet
from refet import qaqc


s_args = {'elev': 1208.5, 'lat': 39.4575, 'zw': 3.0}


def make_daily(**overrides):
    args = dict(
        tmin=19.25, tmax=39.33, tdew=9.91, rs=28.22, uz=2.15,
        zw=s_args['zw'], elev=s_args['elev'], lat=s_args['lat'], doy=182)
    args.update(overrides)
    return refet.Daily(**args)


def test_check_clean_day():
    flags = qaqc.check(make_daily())
    assert set(flags) == {
        'tmax_not_above_tmin', 'ea_above_saturation', 'rs_negative',
        'rs_above_rso', 'ea_negative', 'uz_negative', 'uz_above_max'}
    assert not any(f.any() for f in flags.values())


def test_check_flags_inverted_temps():
    flags = qaqc.check(make_daily(tmin=39.33, tmax=19.25, tdew=5.0))
    assert flags['tmax_not_above_tmin'].all()


def test_check_flags_supersaturation():
    # Dew point far above tmax forces ea > es(tmax)
    flags = qaqc.check(make_daily(tdew=45.0))
    assert flags['ea_above_saturation'].all()


def test_check_flags_rs_above_rso():
    flags = qaqc.check(make_daily(rs=60.0))
    assert flags['rs_above_rso'].all()


def test_check_flags_negative_and_extreme():
    flags = qaqc.check(make_daily(rs=-1.0, uz=45.0))
    assert flags['rs_negative'].all()
    assert flags['uz_above_max'].all()
    flags = qaqc.check(make_daily(uz=-0.1))
    assert flags['uz_negative'].all()


def test_flag_functions_shapes():
    rs = np.array([10.0, 20.0, 30.0])
    rso = np.array([25.0, 25.0, 25.0])
    flag = qaqc.flag_rs_above_rso(rs, rso)
    assert flag.tolist() == [False, False, True]
    assert qaqc.flag_range(np.array([-1.0, 50.0, 200.0]), 0, 100).tolist() == \
        [True, False, True]


def test_counts():
    flags = {'a': np.array([True, False, True]), 'b': np.array([False])}
    assert qaqc.counts(flags) == {'a': 2, 'b': 0}


def test_check_hourly_matches_vpd_clamp():
    """The hourly saturation flag marks exactly the clamped-VPD time steps"""
    pd = pytest.importorskip('pandas')
    import os
    data = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    df = pd.read_csv(os.path.join(data, 'FALN_Agrimet_hourly_raw_2015.csv'))
    h = refet.Hourly(
        tmean=df['OB'].values, tdew=df['TP'].values, rs=df['SI'].values,
        uz=df['WS'].values, zw=3.0, elev=1208.5, lat=39.4575,
        lon=-118.77388, doy=df['DAY'].values * 0 + 182,
        time=df['HOUR'].values,
        input_units={'tmean': 'F', 'tdew': 'F', 'rs': 'Langleys', 'uz': 'mph'})
    flags = qaqc.check(h)
    # 164 hours in the 2015 Fallon record have tdew above tmean
    assert qaqc.counts(flags)['ea_above_saturation'] == 164

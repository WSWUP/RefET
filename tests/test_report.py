import datetime as dt
import os

import pytest

pd = pytest.importorskip('pandas')

from refet.report import ConfigError, build_report
from refet.report.facts import month_day, month_day_year
from refet.report.render import apply_sections, make_ticks, render

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


# ------------------------------------------------------------ date formatting
def test_month_day_windows_safe():
    # strftime('%-d') raises on native Windows; these helpers must not
    assert month_day(dt.date(2015, 7, 1)) == 'Jul 1'
    assert month_day(dt.date(2015, 12, 25)) == 'Dec 25'
    assert month_day_year(dt.date(2015, 7, 1)) == 'Jul 1, 2015'


def test_make_ticks_short_range():
    dates = [f'2015-07-{d:02d}' for d in range(1, 6)]
    ticks = make_ticks(dates)
    assert ticks[0] == {'i': 0, 'label': 'Jul 1'}


def test_make_ticks_full_year():
    start = dt.date(2015, 1, 1)
    dates = [(start + dt.timedelta(days=i)).isoformat() for i in range(365)]
    ticks = make_ticks(dates)
    assert ticks[0]['label'] == 'Jan'
    assert len(ticks) == 12


# ------------------------------------------------------------ render contract
def test_render_substitutes():
    out = render('<p>{{a}}</p><script>const D=__DATA__;</script>',
                 {'a': 'hello'}, '{"x":1}')
    assert '<p>hello</p>' in out
    assert 'const D={"x":1};' in out


def test_render_undefined_placeholder_fails():
    with pytest.raises(ConfigError, match='undefined placeholders'):
        render('{{nope}} __DATA__', {}, '{}')


def test_render_unused_fact_fails():
    with pytest.raises(ConfigError, match='never used'):
        render('{{a}} __DATA__', {'a': '1', 'orphan': '2'}, '{}')


def test_render_js_reference_counts_as_used():
    out = render('{{a}} <script>let x = F_.b; const D=__DATA__;</script>',
                 {'a': '1', 'b': '2'}, '{}')
    assert 'F_.b' in out


def test_render_missing_data_placeholder_fails():
    with pytest.raises(ConfigError, match='__DATA__'):
        render('{{a}}', {'a': '1'}, '{}')


def test_apply_sections_gating():
    tpl = 'always<!--if maybe--> sometimes<!--endif-->'
    assert apply_sections(tpl, {'maybe': True}) == 'always sometimes'
    assert apply_sections(tpl, {'maybe': False}) == 'always'


def test_apply_sections_unknown_flag_fails():
    with pytest.raises(ConfigError, match='no such flag'):
        apply_sections('<!--if zz-->x<!--endif-->', {})


def test_render_gated_fact_still_counts_as_used():
    raw = '{{a}}<!--if v-->{{b}}<!--endif--> __DATA__'
    tpl = apply_sections(raw, {'v': False})
    out = render(tpl, {'a': '1', 'b': '2'}, '{}', usage_source=raw)
    assert out == '1 {}'


# --------------------------------------------------- offline end-to-end build
def faln_2015_config(tmp_path, with_reference=True):
    data = DATA.replace(os.sep, '/')
    ref_daily = (f"reference = '{data}/FALN_Agrimet_daily_raw_2015.out'"
                 if with_reference else '')
    ref_hourly = (f"reference = '{data}/FALN_Agrimet_hourly_raw_2015.out'"
                  if with_reference else '')
    station = '' if with_reference else (
        '[station]\nlat = 39.4575\nlon = -118.77388\nelev = 1208.5\nzw = 3.0\n')
    text = f"""
title = "Fallon, Nevada"
name = "Fallon, NV (FALN) AgriMet"
network = "AgriMet"
start = 2015-01-01
end = 2015-12-31
timezone = "US/Pacific"
method = "refet"
rso_type = "full"

[source]
type = "csv"

{station}
[daily]
csv = '{data}/FALN_Agrimet_daily_raw_2015.csv'
{ref_daily}

[daily.columns]
tmin = "MN"
tmax = "MX"
tdew = "YM"
uz = "UA"
rs = "SR"

[hourly]
csv = '{data}/FALN_Agrimet_hourly_raw_2015.csv'
{ref_hourly}

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
"""
    path = tmp_path / 'faln_2015.toml'
    path.write_text(text, encoding='utf-8')
    return str(path)


def test_build_report_faln_2015(tmp_path):
    """Full offline build of the validated Fallon 2015 report."""
    config = faln_2015_config(tmp_path)
    output = str(tmp_path / 'report.html')
    payload, facts, out_path = build_report(config, output_path=output)

    html = open(output, encoding='utf-8').read()
    assert 'Fallon' in html
    assert '__DATA__' not in html
    assert '{{' not in html

    # The report's own validation panel is the regression gate: computed
    # daily ETr must track the Ref-ET .out values
    v = payload['validation']
    assert v['daily_etr']['n'] >= 360
    assert v['daily_etr']['rmse'] < 0.01
    assert v['daily_eto']['max_abs'] < 0.05
    assert v['hourly_etr']['rmse'] < 0.05
    assert payload['flags']['validation'] is True

    # The 2015 record has one day with missing inputs
    assert facts['n_days'] == '364'
    assert facts['period_phrase'] == 'A full year'
    assert len(payload['heatmap']) == 24
    assert len(payload['dates']) == 365


def test_build_report_without_reference(tmp_path):
    """No reference files: validation panel omitted, build still succeeds."""
    config = faln_2015_config(tmp_path, with_reference=False)
    output = str(tmp_path / 'report.html')
    payload, facts, out_path = build_report(config, output_path=output)
    assert payload['validation'] is None
    assert payload['flags']['validation'] is False
    assert 'val_rmse' not in facts
    assert os.path.exists(output)


def test_build_report_bad_config(tmp_path):
    path = tmp_path / 'bad.toml'
    path.write_text('title = "x"\n', encoding='utf-8')
    with pytest.raises(ConfigError, match='missing required key'):
        build_report(str(path))

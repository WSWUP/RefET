"""Build an interactive HTML report of reference ET for a station.

Reads a TOML station config, computes daily and hourly ASCE reference ET
with this library, optionally validates against Ref-ET .out files, and
renders the result into a single self-contained HTML page::

    refet-report stations/faln_2015.toml

See the example configs in the repository's examples/stations/ directory.
"""
import argparse
import json
import os
import sys
from importlib import resources


def default_template():
    return (resources.files('refet.report') / 'template.html'
            ).read_text(encoding='utf-8')


def build_report(config_path, template_path=None, output_path=None,
                 json_path=None):
    """Build one report; returns (payload, facts, output_path)."""
    from .compute import build_payload, load_sources
    from .config import load_config
    from .render import apply_sections, render

    cfg = load_config(config_path)
    daily_raw, hourly_raw, source_note = load_sources(cfg)
    payload, facts = build_payload(cfg, daily_raw, hourly_raw, source_note)

    if template_path:
        with open(template_path, encoding='utf-8') as f:
            raw_template = f.read()
    else:
        raw_template = default_template()

    output = output_path or (
        os.path.splitext(os.path.basename(config_path))[0] + '_report.html')
    data_json = json.dumps(payload)
    template = apply_sections(raw_template, payload['flags'])
    html = render(template, facts, data_json, usage_source=raw_template)
    with open(output, 'w', encoding='utf-8') as f:
        f.write(html)
    if json_path:
        with open(json_path, 'w', encoding='utf-8') as f:
            f.write(data_json)
    return payload, facts, output


def main(argv=None):
    p = argparse.ArgumentParser(
        prog='refet-report',
        description='Build a self-contained HTML reference ET report from a '
                    'TOML station config.')
    p.add_argument('config', help='TOML station config')
    p.add_argument('--template', help='override the built-in report template')
    p.add_argument('--output', help='output HTML path '
                   '(default: <config name>_report.html in the current dir)')
    p.add_argument('--json', help='also write the raw data payload here')
    args = p.parse_args(argv)

    try:
        import pandas  # noqa: F401
    except ImportError:
        sys.exit('error: the report builder needs pandas; '
                 'install with: pip install refet[report]')

    from .config import ConfigError

    try:
        payload, facts, output = build_report(
            args.config, args.template, args.output, args.json)
    except ConfigError as e:
        sys.exit(f'error: {e}')

    print(f"{facts['station_name']}, {facts['start_date']} to "
          f"{facts['end_date']}: ETr {facts['etr_total']} mm / "
          f"ETo {facts['eto_total']} mm over {facts['n_days']} days and "
          f"{facts['n_hours']} hours")
    if payload['validation']:
        for key in ('daily_etr', 'daily_eto', 'hourly_etr'):
            s = payload['validation'].get(key)
            if s:
                print(f"  {key:<11} n={s['n']:<6} rmse={s['rmse']:.4f} "
                      f"max={s['max_abs']:.4f} bias={s['bias']:+.4f} mm")
    else:
        print('  no reference files configured; validation panel omitted')
    print(f'Wrote {output}')


if __name__ == '__main__':
    main()

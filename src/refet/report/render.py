"""Strict template rendering.

Every number in the page's prose is a derived fact. The template writes them
as ``{{placeholder}}``; page scripts read them as ``F_.name``. Both
directions are checked: a placeholder with no fact means the prose states
something nothing derives; a fact nothing references means the narrative
silently dropped a number it used to report. Either is a build failure.
"""
import datetime as dt
import re

from .config import ConfigError
from .facts import month_day

PLACEHOLDER = re.compile(r'\{\{\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\}\}')
# Facts the page's own scripts read at runtime, e.g. F_.peak_date
JS_REFERENCE = re.compile(r'\bF_\.([a-zA-Z_][a-zA-Z0-9_]*)\b')
# Optional template region: <!--if validation--> ... <!--endif-->
SECTION = re.compile(r'[ \t]*<!--\s*if\s+(\w+)\s*-->\n?(.*?)'
                     r'[ \t]*<!--\s*endif\s*-->\n?', re.S)


def make_ticks(dates):
    """Axis tick positions and labels, chosen from the span of the record."""
    n = len(dates)
    parsed = [dt.date.fromisoformat(d) for d in dates]
    multi_year = parsed[0].year != parsed[-1].year
    ticks = []
    if n <= 16:
        step = 1 if n <= 8 else 2
        for i in range(0, n, step):
            ticks.append({'i': i, 'label': month_day(parsed[i])})
    elif n <= 70:
        for i in range(0, n, 7):
            ticks.append({'i': i, 'label': month_day(parsed[i])})
    else:
        # Month starts, thinned so labels never crowd
        starts = [i for i, d in enumerate(parsed) if d.day == 1] or [0]
        every = max(1, round(len(starts) / 12))
        for i in starts[::every]:
            d = parsed[i]
            label = d.strftime('%b')
            if multi_year and (d.month == 1 or i == starts[0]):
                label += d.strftime(' %Y')
            ticks.append({'i': i, 'label': label})
    return ticks


def apply_sections(template, flags):
    """Keep or drop <!--if name--> regions before placeholders are resolved.

    Without this, a report with no reference data would fail the placeholder
    check on validation numbers that legitimately do not exist.
    """
    def sub(m):
        name = m.group(1)
        if name not in flags:
            raise ConfigError(f'template has <!--if {name}--> but no such flag')
        return m.group(2) if flags[name] else ''
    return SECTION.sub(sub, template)


def render(template, facts, data_json, usage_source=None):
    """Substitute {{facts}} and the data payload, strictly.

    Usage is measured against the template *before* optional sections are
    dropped (``usage_source``), so a fact referenced only inside a gated-off
    region still counts as used rather than looking like dead weight.
    """
    scan = usage_source if usage_source is not None else template
    used = set(PLACEHOLDER.findall(scan)) | set(JS_REFERENCE.findall(scan))
    unknown = sorted(set(PLACEHOLDER.findall(template)) - set(facts))
    if unknown:
        raise ConfigError(
            f'template uses undefined placeholders: {unknown}\n'
            f'Add them to derive_facts() or remove them from the template.')
    unused = sorted(set(facts) - used)
    if unused:
        raise ConfigError(
            f'derived facts never used by the template: {unused}\n'
            f'Reference them as {{{{name}}}} in the prose or F_.name in the '
            f'page script, or drop them from derive_facts().')
    out = PLACEHOLDER.sub(lambda m: facts[m.group(1)], template)
    if '__DATA__' not in out:
        raise ConfigError('template has no __DATA__ placeholder for the payload')
    return out.replace('__DATA__', data_json)

"""Self-contained HTML reference ET reports (``pip install refet[report]``).

The public entry points are :func:`build_report` and the ``refet-report``
console script. See examples/stations/ in the repository for worked configs.
"""
from .cli import build_report
from .config import ConfigError

__all__ = ['build_report', 'ConfigError']

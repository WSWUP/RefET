from .daily import Daily
from .hourly import Hourly

try:
    from importlib import metadata
except ImportError:  # for Python<3.8
    import importlib_metadata as metadata

__version__ = metadata.version("refet")

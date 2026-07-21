from importlib import metadata

from . import calcs
from . import estimate
from . import qaqc
from . import units
from .daily import Daily
from .hourly import Hourly

__version__ = metadata.version("refet")

__all__ = ['Daily', 'Hourly', 'calcs', 'estimate', 'qaqc', 'units', '__version__']

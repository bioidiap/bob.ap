from ._library import *
from . import version
from .version import module as __version__

# gets sphinx autodoc done right - don't remove it
__all__ = [_ for _ in dir() if not _.startswith('_')]

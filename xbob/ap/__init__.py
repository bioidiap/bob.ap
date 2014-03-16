from ._library import __version__
from ._library import *

# gets sphinx autodoc done right - don't remove it
__all__ = [k for k in dir() if not k.startswith('_')]
del k

import logging
import warnings

from .HuMPI import *

try:
	logging.lastResort
except AttributeError:
	logging.getLogger(__name__).addHandler(logging.StreamHandler())



"""
Python module containing calibration data and lookup functions for standard
thermocouples of types B, C, D, E, G, J, K, N, R, S, T.

You can access the lookup table objects like so::

    import thermocouples_reference
    typeK = thermocouples_reference.thermocouples['K']

Find the latest version at https://pypi.python.org/pypi/thermocouples_reference

===DISCLAIMER===

This module is provided for educational purposes. For any real-world
process, I strongly recommend that you check the output of this module
against a known good standard.

I make no warranties as to the accuracy of this module, and shall
not be liable for any damage that may result from errors or omissions.

================

"""

__author__    = "User:Nanite @ wikipedia"
__copyright__ = "public domain"
__version__   = "0.20"

from . import units
from . import function_types
from . import source_NIST
from . import source_ASTM
from . import source_OMEGA

# assemble thermocouples list.
# note the order: type G from OMEGA source is replaced by the ASTM one.
thermocouples = dict()
thermocouples.update(source_OMEGA.thermocouples)
thermocouples.update(source_NIST .thermocouples)
thermocouples.update(source_ASTM .thermocouples)

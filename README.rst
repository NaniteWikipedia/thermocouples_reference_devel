=======================
thermocouples_reference
=======================

Python module containing calibration data and lookup functions for standard
`thermocouples`_ of types **B**, **C**, **D**, **E**, **G**, **J**, **K**,
**M**, **N**, **P**, **R**, **S**, **T**, and some less standard types too.

Usage and examples
------------------

Below, the first computation shows that the type K thermocouple
emf at 42 °C, with reference junction at 0 °C, is 1.694 mV
(`compare to NIST table`_); the second calculation shows how passing
in an array applies the function for each element, in the style of numpy:

  >>> from thermocouples_reference import thermocouples
  >>> typeK = thermocouples['K']
  >>> typeK
  <Type K thermocouple reference (-270.0 to 1372.0 degC)>
  >>> typeK.emf_mVC(42, Tref=0)
  1.6938477049901346
  >>> typeK.emf_mVC([-3.14159, 42, 54], Tref=0)
  array([-0.12369326,  1.6938477 ,  2.18822176])

An inverse lookup function is provided that you can use to get a temperature
out of a measured voltage, including cold junction compensation effects.
If we put our type K thermocouple into a piece of spam and we read 1.1 mV,
using our voltmeter at room temperature (23 °C), then the spam is at
50 °C. [1]_

  >>> typeK.inverse_CmV(1.1, Tref=23.0)
  49.907928030075773
  >>> typeK.emf_mVC(49.907928030075773, Tref=23.0) # check result
  1.1000000000000001

The functions are called ``emf_mVC`` and ``inverse_CmV`` just to remind you
about the units of voltage and temperature. Other temperature units are
supported as well:

==================   ==========   ==============
 Temperature unit    EMF lookup   Inverse lookup
==================   ==========   ==============
degrees Celsius      .emf_mVC     .inverse_CmV
degrees Fahrenheit   .emf_mVF     .inverse_FmV
kelvins              .emf_mVK     .inverse_KmV
degrees Rankine      .emf_mVR     .inverse_RmV
==================   ==========   ==============

You can also compute derivatives of the emf function. These are functional
derivatives, not finite differences. The Seebeck coefficients of chromel
and alumel differ by 42.00 μV/°C, at 687 °C:

  >>> typeK.emf_mVC(687,derivative=1)
  0.041998175982382979

Data sources
------------

Readers may be familiar with thermocouple lookup tables (`example table`_).
Such tables are computed from standard reference functions, generally
piecewise polynomials. [2]_ This module contains the source polynomials
*directly*, and so in principle it is more accurate than any lookup table.
Lookup tables also often also include approximate polynomials for temperature
lookup based on a given compensated emf value. Such inverse polynomials are
*not included* in this module; rather, the inverse lookup is based on
numerically searching for a solution on the exact emf function.

For any thermocouple object, information about calibration and source is
available in the repr() of the .func attribute:

    >>> typeK.func
    <piecewise polynomial+gaussian, domain -270.0 to 1372.0 in degC, output in mV; 
    ITS-90 calibrated, from NIST SRD 60, type K>

The data sources are:

- Types B, E, J, K, N, R, S, T
  use coefficients from `NIST's website`_, and are calibrations
  to the `ITS-90`_ scale. [3]_
- Types G, M, P, and non-lettered types Au-Pt, Au-Pd, AuFe 0.07,
  IrRh 40-0, PtMo 5-0.1, PtRh 40-20
  use coefficients from `ASTM E 1751-00`_ and are calibrations to ITS-90.
- Types C, D [4]_
  use coefficients found from a publication of OMEGA Engineering
  Inc., and are calibrations to `IPTS-68`_ scale. [5]_

Graphs of functions (if you don't see anything, see
`low temperature types here`_, `intermediate temperature types here`_, and
`high temperature types here`_):

.. image:: https://upload.wikimedia.org/wikipedia/commons/f/f8/Low_temperature_thermocouples_reference_functions.svg
.. image:: https://upload.wikimedia.org/wikipedia/commons/9/95/Intermediate_temperature_thermocouples_reference_functions.svg
.. image:: https://upload.wikimedia.org/wikipedia/commons/c/c3/High_temperature_thermocouples_reference_functions.svg

Requirements
------------

- ``numpy``
- ``scipy`` (optional, only needed for inverse lookup)
- ``python2`` or ``python3`` languages

Installation
------------

Recommended installation is via pip. First, `install pip`_. Then::

    pip install thermocouples_reference --user

(Remove the ``--user`` option if you are superuser and want to install
system-wide.)

Disclaimer
----------
This module is provided for educational purposes. For any real-world
process, I strongly recommend that you check the output of this module
against a known good standard.

I make no warranties as to the accuracy of this module, and shall
not be liable for any damage that may result from errors or omissions.

.. _thermocouples: https://en.wikipedia.org/wiki/Thermocouple
.. _emf reference function: https://en.wikipedia.org/wiki/Thermocouple#Thermocouple_characteristic_function
.. _install pip: http://www.pip-installer.org/en/latest/installing.html
.. _compare to NIST table: http://srdata.nist.gov/its90/download/type_k.tab
.. _low temperature types here: http://commons.wikimedia.org/wiki/File:Low_temperature_thermocouples_reference_functions.svg
.. _intermediate temperature types here: http://commons.wikimedia.org/wiki/File:Intermediate_temperature_thermocouples_reference_functions.svg
.. _high temperature types here: http://commons.wikimedia.org/wiki/File:High_temperature_thermocouples_reference_functions.svg
.. _NIST's website: http://srdata.nist.gov/its90/main/
.. _example table: http://srdata.nist.gov/its90/download/type_k.tab
.. _ITS-90: https://en.wikipedia.org/wiki/International_Temperature_Scale_of_1990
.. _ASTM E 1751-00: http://www.google.com/search?q=ASTM+E1751
.. _IPTS-68: http://www.bipm.org/en/si/history-si/temp_scales/ipts-68.html
.. [1] This is the optimal temperature for spam. Always make sure your
       spam reads around 1.1 millivolt and you'll have a tasty treat.
.. [2] A notable exception is NIST's type K curve which uses a polynomial plus
       gaussian. The gaussian conveniently captures a wiggle in the Seebeck
       coefficient of alumel, that happens around 130 °C.
.. [3] The ITS-90 value *T*\ :sub:`90` is believed to track the true
       thermodynamic temperature *T* very closely. 
       The error *T* − *T*\ :sub:`90` is quite small, of order 0.01 K for
       everyday conditions (up to about 200 °C), rising to around 0.05 K up
       at 1000 °C, and increasing even further after that. See
       `Supplementary Information for the ITS-90`_. Generally your
       thermocouple accuracy will be more limited by manufacturing variations
       and by degradation of the metals in the thermal gradient region.
.. [4] An extra type G IPTS68 curve from the same source is available in
       ``thermocouples_reference.source_OMEGA.thermocouples``. The type G in
       the main ``thermocouples_reference.thermocouples`` contains the ASTM
       curve which is ITS-90 calibrated.
.. [5] IPTS-68 reads higher than ITS-90 by about 1 °C at high temperatures
       around 2000 °C. See `Supplementary Information for the ITS-90`_
       (specifically Fig. 5 in the Introduction) for more information about
       the difference.

.. _Supplementary Information for the ITS-90: http://www.bipm.org/en/publications/mep_kelvin/its-90_supplementary.html
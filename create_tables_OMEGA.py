#!/usr/bin/python
"""
This module creates python code for thermocouple tables.

This module uses thermocouple reference functions,
extracted from the coefficient tables downloaded from
    "Tungsten-Rhenium Thermocouples Calibration Equivalents"
    http://www.omega.com/temperature/z/pdf/z202.pdf
Note: these are probably IPTS68 calibrated

Disclaimers:
(Author) I make no warranties as to the accuracy of this module, and shall
        not be liable for any damage that may result from errors or omissions.
"""

__author__    = "User:Nanite @ wikipedia"
__copyright__ = "public domain"

import numpy as np
import scipy.linalg

def polyshift_mtx(deg, mult, shift):
    """
    Compute matrix for polynomial shift and multiply.
    This is a helper for mapping Fahrenheit to Celsius polynomials.

    deg : int
        Polynomial degree.
    mult : float
        x-multiplier value (see below for meaning)
    shift : float
        x-shift value (see below for meaning)
    
    The output is a matrix used to transform polynomial coefficents.
    Note that it is in ascending power, so for use with numpy.polyval() you need
    to reverse it as below.
    The shift and multiply are defined so that for any coefs,x,shift,mult
        deg = len(coefs)-1
        mtx = polyshift_mtx(deg, shift, mult)[::-1,::-1]
        coefs2 = mtx.dot(coefs)
        polyval(coefs2, x) == polyval(coefs, mult*(x+shift))
    """
    deg = int(deg)
    shift = float(shift)
    mult = float(mult)

    multpowers = np.arange(deg+1, dtype=float)
    multpowers[1:] = mult ** multpowers[1:]
    multpowers[0] = 1.

    shiftpowers = np.arange(deg+1, dtype=float)
    shiftpowers[1:] = shift ** shiftpowers[1:]
    shiftpowers[0] = 1.
    
    mtx = scipy.linalg.pascal(deg+1, kind='upper', exact=False)
    for col in range(deg+1):
        mtx[col::-1,col] *= multpowers[col]*shiftpowers[:col+1]
    
    return mtx

_conv_CF = polyshift_mtx(5,9./5.,32.*5./9.)[::-1,::-1]

typeG = _conv_CF.dot(np.array([
        -0.1176051e-16,
         0.2125270e-12,
        -0.1795965e-8,
         0.6783829e-5,
         0.2883146e-3,
        -0.1580014e-1,
    ]))
typeG[-1] = 0. # put zero reference at 0 deg C

typeC = _conv_CF.dot(np.array([
        -0.2616792e-16,
         0.3471851e-12,
        -0.1842722e-8,
         0.3956443e-5,
         0.7190027e-2,
        -0.234471,
    ]))
typeC[-1] = 0. # put zero reference at 0 deg C

output = '''\
"""    
This module contains thermocouple reference functions for types C,D,G.

You can access the lookup table objects like so:
    typeC = <this module>.thermocouples['C']

Thermocouple reference functions are polynomials taken from
    "Tungsten-Rhenium Thermocouples Calibration Equivalents"
    http://www.omega.com/temperature/z/pdf/z202.pdf

Note on calibration
-------------------
Curves C, G are almost certainly calibrations to IPTS-68,
as suggested by the source PDF.
Observe that we can compare the type G curve here to
the type G from ASTM functions:
    from pylab import *
    from thermocouples_reference import *
    G90 = source_ASTM.thermocouples['G']
    T90 = linspace(0,2310,1001)
    G68 = source_OMEGA.thermocouples['G']
    T68 = [G68.inverse_CmV(e) for e in G90.emf_mVC(T90)]
    figure() ; plot(T90+273.15,T90-T68)
    xlim(0,4000) ; ylim(-2.5,0.5)
    xlabel('$T_{90}$ (K)')
    ylabel('$T_{90} - T_{68}$ (K)')
    show()
The resulting graph appears similar to the published
ITS-90 vs. IPTS-68 difference, see 
http://www.bipm.org/en/publications/mep_kelvin/its-90_supplementary.html
    Figure 5 of Introduction

The difference between IPTS-68 and ITS-90 scales is at worst
* 0.4 deg C over the range 0 to 1400 deg C;
* At 1500 deg C, IPTS-68 reads 0.44 deg C higher
* At 2000 deg C, IPTS-68 reads 0.72 deg C higher
* At 2500 deg C, IPTS-68 reads 1.07 deg C higher
Anyway at the manufacturing variations in WRe thermocouples
are somewhere around +/-4 to +/-20 deg C.

Curve D is also very probably IPTS68.

Disclaimers
-----------
(Author) I make no warranties as to the accuracy of this module, and shall
        not be liable for any damage that may result from errors or omissions.

(Note: This module is generated code from create_tables_OMEGA.py.)
"""

__copyright__ = "public domain"

import numpy as np
from .function_types import Thermocouple_Reference, Polynomial_Gaussian_Piecewise_Function

_source = 'OMEGA Inc. z202.pdf'
_cal = 'IPTS-68'

thermocouples = {
# The coefficients of this polynomial have been converted from
# a Fahrenheit polynomial with 7 significant figures, and have
# been stored here with longer precision to avoid further inaccuracy.
'G':Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function(
    [[0.,2315.,
np.'''+np.array_repr(typeG,precision=16)+''',
    None]],'C','mV', calibration=_cal, source=_source+', type G'),
    ttype='Type G',
    composition='W - 74W,26Re'),

# The coefficients of this polynomial have been converted from
# a Fahrenheit polynomial with 7 significant figures, and have
# been stored here with longer precision to avoid further inaccuracy.
'C':Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function(
    [[0.,2315.,
np.'''+np.array_repr(typeC,precision=16)+''',
    None]],'C','mV', calibration=_cal, source=_source+', type C'),
    ttype='Type C',
    composition='95W,5Re - 74W,26Re'),

'D':Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function(
    [[0.,783., np.array([
        -1.4240735e-15,
         7.9498033e-12,
        -1.8464573e-8,
         2.0592621e-5,
         9.5685256e-3,
         0.,
    ]), None],
    [783.,2320., np.array([
        -7.9026726e-16,
         5.3743821e-12,
        -1.4935266e-8,
         1.8666488e-5,
         9.9109462e-3,
         0.,
    ]), None],
    ],'C','mV', calibration=_cal, source=_source+', type D'),
    ttype='Type D',
    composition='97W,3Re - 75W,25Re'),
}

del _source
del _cal

#end of module
'''

print(output)



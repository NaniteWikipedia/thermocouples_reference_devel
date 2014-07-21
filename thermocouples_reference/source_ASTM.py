"""
This module contains thermocouple reference functions from ASTM E 1751-00,
according to which none of these have assigned letters.

In fact, a few of them do have "unofficial" assigned letters. The website
http://www.capgo.com/Resources/Temperature/Thermocouple/Thermocouple.html
(as of 2013 December) has a nice table assigning some letters to unusual
types contained herein.

You can access the lookup table objects like so:
    typeM = <this module>.thermocouples['M']

About these thermocouples
-------------------------
What are these strange thermocouples used for? Some useful information
about their purposes can be found in the following literature

* Thermoelectricity: Theory, Thermometry, Tool, Issue 852 by Daniel D. Pollock
* fluke.com : search for 5629 Gold Platinum Thermocouple

A summary of properties:

==========  ========  ======
Identifier   Table    Notes
==========  ========  ======
G           Table  1  W-Re system, very high temp.
P           Table  3  Platinel - noble metal mix that mimics type K emf.
AuFe 0.07   Table  5  Fe-doped Au vs chromel, used in cryogenics.
PtMo 5-01   Table  7  Pt-Mo system, used in nuclear reactors.
                      Low drift from transmutation.
PtRh 40-20  Table  9  Pt-Rh system (like types B,R,S) but goes to 1800C
M           Table 11  Ni-Mo/Co system, similar to other Ni-based (e.g. K)
                      but limit ~1500C is a bit better.
IrRh 40-0   Table 13  Ir-Rh system, goes up to ~2000C in inert.
Au-Pt       Table 15  Pure metals: amazing stability, but limited
                      because gold melts at 1064C.
Pt-Pd       Table 17  Also pure. Palladium melts at 1555C.
==========  ========  ======

Disclaimers
-----------
(Author) I make no warranties as to the accuracy of this module, and shall
        not be liable for any damage that may result from errors or omissions.
"""

__copyright__ = "public domain"

import numpy as np
from .function_types import Thermocouple_Reference, Polynomial_Gaussian_Piecewise_Function

_source = "ASTM E 1751-00"
_cal = "ITS-90"

thermocouples = {
'G': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0., 630.615, np.array([
        -1.7089202e-15,
         4.3850022e-12,
        -1.1393234e-08,
         2.1634754e-05,
         1.2792201e-03,
         0.0000000e+00,
  ]), None],
  [630.615, 2315., np.array([
        -1.5534591e-25,
         1.8120237e-21,
        -8.9888053e-18,
         2.4455012e-14,
        -3.8615222e-11,
         3.1141330e-08,
        -3.6467516e-06,
         9.4962455e-03,
        -1.1064412e+00,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 1"),
    ttype='Type G',
    composition='W - 74W,26Re'),

'P': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0, 746.6, np.array([
        -3.6375467e-15,
         1.4851327e-11,
        -3.4878428e-08,
         3.5175152e-05,
         2.9819716e-02,
         0.0000000e+00,
  ]), None],
  [746.6, 1395, np.array([
        -9.3211269e-18,
         5.4438760e-14,
        -1.2855115e-10,
         1.5424937e-07,
        -1.0570233e-04,
         8.5377200e-02,
        -8.9621838e+00,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 3"),
    ttype='Type P (II)',
    composition='55Pd,31Pt,14Au - 65Au,35Pd'),

'AuFe 0.07': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [-273., 7., np.array([
         6.8263661580e-31,
         1.1010930596e-27,
         7.8225430483e-25,
         3.2146639387e-22,
         8.4287909747e-20,
         1.4636450149e-17,
         1.6829773697e-15,
         1.2272348484e-13,
         4.9063035769e-12,
         4.0432555769e-11,
        -4.5260169888e-09,
        -1.5967928202e-07,
         3.6406179664e-06,
         2.2272367466e-02,
         0.0000000000e+00,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 5"),
    ttype='Chromel-AuFe0.07',
    composition='90Ni,10Cr - Au,0.07(atom%)Fe'),

'PtMo 5-0.1': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0., 491., np.array([
        -2.0186476e-19,
         3.3574252e-16,
        -2.3848950e-13,
         1.0585770e-10,
        -4.3368594e-08,
         2.8410937e-05,
         1.0501456e-02,
         0.0000000e+00,
  ]), None],
  [491., 1600., np.array([
         9.4670862e-24,
        -7.6717268e-20,
         2.6865173e-16,
        -5.3071212e-13,
         6.4615219e-10,
        -4.9920472e-07,
         2.4913353e-04,
        -4.8776479e-02,
         6.8354086e+00,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 7"),
    ttype='PtMo 5-0.1',
    composition='95Pt,5Mo - 99.9Pt,0.1Mo'),

'PtRh 40-20': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0, 951.7, np.array([
        -2.8497160e-22,
         1.0033974e-18,
        -1.5406939e-15,
         1.0382985e-12,
         4.2594137e-10,
         3.9360320e-07,
         3.6246289e-04,
         0.0000000e+00,
  ]), None],
  [951.7, 1888., np.array([
        -1.2619640e-20,
         1.1516280e-16,
        -1.0824710e-12,
         3.6728697e-09,
        -3.9077442e-06,
         3.5246931e-03,
        -9.1201877e-01,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 9"),
    ttype='PtRh 40-20',
    composition='60Pt,40Rh - 80Pt,20Rh'),
# 

'M': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [-50., 370.8, np.array([
        -3.394387900e-19,
        -9.738054601e-17,
         1.846977453e-13,
        -1.025216130e-10,
        -3.142898226e-08,
         4.408522682e-05,
         3.690092195e-02,
         0.000000000e+00,
  ]), None],
  [370.8, 1410., np.array([
         1.027600874e-25,
        -7.864442961e-22,
         2.627522669e-18,
        -5.041679909e-15,
         6.145877457e-12,
        -4.958763813e-09,
         2.650568429e-06,
        -8.846963426e-04,
         2.059913943e-01,
        -1.145582129e+01,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 11"),
    ttype='Type M',
    composition='82Ni,18Mo - 99.2Ni,0.8Co'),

'IrRh 40-0': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0., 630.615, np.array([
        -7.9634082e-23,
         1.5270867e-19,
        -1.0418040e-16,
         2.6762413e-14,
         2.7700591e-12,
        -7.8890504e-09,
         6.9649773e-06,
         3.0870016e-03,
         0.0000000e+00,
  ]), None],
  [630.615, 2110., np.array([
         3.0821886e-20,
        -5.1797037e-16,
         2.7235393e-12,
        -6.0547943e-09,
         5.7455189e-06,
         3.6588615e-03,
        -9.6839082e-02,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 13"),
    ttype='IrRh 40-0',
    composition='60Ir,40Rh - Ir'),

'Au-Pt': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0., 1000., 1e-3*np.array([ # coefs in microvolt
        -2.51672787e-24,
         1.42981590e-20,
        -3.39430259e-17,
         4.56927038e-14,
        -4.24206193e-11,
         3.28711859e-08,
        -2.22998614e-05,
         1.93672974e-02,
         6.03619861e+00,
         0.00000000e+00,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 15"),
    ttype='Au-Pt',
    composition='Au - Pt'),

'Pt-Pd': Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([
  [0., 660.323, 1e-3*np.array([ # coefs in microvolt
        -8.510068e-21,
         2.257823e-17,
        -1.268514e-14,
        -2.012523e-11,
         2.992243e-08,
        -9.602271e-06,
         4.610494e-03,
         5.296958e+00,
         0.000000e+00,
  ]), None],
  [660.323, 1500., 1e-3*np.array([
        -1.3570737e-15,
         9.5627366e-12,
        -2.6901509e-08,
         3.6361700e-05,
        -1.5793515e-02,
         1.0182545e+01,
        -4.9771370e+02,
  ]), None],
  ],'C','mV',calibration=_cal,source=_source+", Table 17"),
    ttype='Pt-Pd',
    composition='Pt - Pd'),

}

del _source
del _cal

#end of module
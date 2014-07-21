#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
"""
This module creates python code for thermocouple tables.

This module uses thermocouple reference functions,
taken from the coefficient tables downloaded from
    http://srdata.nist.gov/its90/download/allcoeff.tab
the NIST Standard Reference Database 60, Version 2.0 (Web Version)
    (see http://srdata.nist.gov/its90/main/ for more information)

The output of this file, but not this file itself, is distributed with the package.

Disclaimers:
(Author) I make no warranties as to the accuracy of this module, and shall
        not be liable for any damage that may result from errors or omissions.
"""

__author__    = "User:Nanite @ wikipedia"
__copyright__ = "public domain"

import numpy as np

source = "NIST SRD 60"
cal = "ITS-90"

def Thermocouple_from_NIST_text(text):
    """ Helper function for turning NIST table text into Thermocouple object. """
    lines = text.splitlines()
    ttype = '?'
    family = 'None'
    table = []
    i=0
    while(i<len(lines)):
        l = lines[i].lstrip().split()
        if l[0] == "name:":
            family = "'"+lines[i].lstrip()[6:]+"'"
        elif l[0] == "type:":
            ttype = l[1]
#        elif l[0] == "temperature" and l[1] == "units:":
#            Tunits = l[2]
#        elif l[0] == "emf" and l[1] == "units:":
#            emfunits = l[2]
        elif l[0] == "range:":
            num = int(l[3])
            entry = [l[1].rstrip(','),
                     l[2].rstrip(','),
                     lines[i+num+1:i:-1], 'None']
            table.append(entry)
            i += num
        elif l[0] == "exponential:":
            if ttype != 'K':
                import warnings
                warnings.warn("Thermocouple_from_NIST_text loading warning: exponential entry encountered for non-K thermocouple type. May be misinterpreted.")
            a0str, a0 = lines[i+1].split('=')
            a1str, a1 = lines[i+2].split('=')
            a2str, a2 = lines[i+3].split('=')
            if a0str.strip() != 'a0' or a1str.strip() != 'a1' or a2str.strip() != 'a2':
                raise ValueError("Invalid exponential entry")
            a0 = a0.strip().lower()
            a1 = a1.strip().lower()
            a2 = a2.strip().lower()
            table[-1][-1] = '['+a0+','+a1+','+a2+']'
            i += 3
        i+=1
    
    outstring = "Thermocouple_Reference(Polynomial_Gaussian_Piecewise_Function([\n"
    for tmin, tmax, numlist, expo in table:
        outstring += "  ["+tmin+", "+tmax+", np.array([\n"
        for num in numlist:
            num = num.strip().lower()
            if num[0] != '-':
                num = ' '+num
            outstring += (" "*8)+num+",\n"
        outstring += "  ]), "+expo+"],\n"
    outstring += "  ],'C','mV',\n"
    outstring += "       calibration='{}',\n".format(cal)
    outstring += "       source='{}, type {}'),\n".format(source,ttype)
    outstring += "     ttype='Type {}',\n".format(ttype)
    outstring += "     composition='{}')".format(comps[ttype])

    return ttype, outstring

comps = {
    'B': '70Pt,30Rh - 94Pt,6Rh',
    'R': '87Pt,13Rh - Pt',
    'S': '90Pt,10Rh - Pt',
    'E': '90Ni,10Cr - 55Cu,45Ni',
    'J': 'Fe - 55Cu,45Ni',
    'K': '90Ni,10Cr - 95Ni,2Mg,2Al,1Si',
    'N': '84.1Ni,14.4Cr,1.4Si,0.1Mg - 95.6Ni,4.4Si',
    'T': 'Cu - 55Cu,45Ni',
    }

nist_texts = ["""\
name: reference function on ITS-90
type: B
temperature units: °C
emf units: mV
range: 0.000, 630.615, 6
  0.000000000000E+00
 -0.246508183460E-03
  0.590404211710E-05
 -0.132579316360E-08
  0.156682919010E-11
 -0.169445292400E-14
  0.629903470940E-18
range: 630.615, 1820.000, 8
 -0.389381686210E+01
  0.285717474700E-01
 -0.848851047850E-04
  0.157852801640E-06
 -0.168353448640E-09
  0.111097940130E-12
 -0.445154310330E-16
  0.989756408210E-20
 -0.937913302890E-24
""","""\
name: reference function on ITS-90
type: E
temperature units: °C
emf units: mV
range: -270.000, 0.000, 13
  0.000000000000E+00
  0.586655087080E-01
  0.454109771240E-04
 -0.779980486860E-06
 -0.258001608430E-07
 -0.594525830570E-09
 -0.932140586670E-11
 -0.102876055340E-12
 -0.803701236210E-15
 -0.439794973910E-17
 -0.164147763550E-19
 -0.396736195160E-22
 -0.558273287210E-25
 -0.346578420130E-28
range: 0.000, 1000.000, 10
  0.000000000000E+00
  0.586655087100E-01
  0.450322755820E-04
  0.289084072120E-07
 -0.330568966520E-09
  0.650244032700E-12
 -0.191974955040E-15
 -0.125366004970E-17
  0.214892175690E-20
 -0.143880417820E-23
  0.359608994810E-27
""","""\
name: reference function on ITS-90
type: J
temperature units: °C
emf units: mV
range: -210.000, 760.000, 8
  0.000000000000E+00
  0.503811878150E-01
  0.304758369300E-04
 -0.856810657200E-07
  0.132281952950E-09
 -0.170529583370E-12
  0.209480906970E-15
 -0.125383953360E-18
  0.156317256970E-22
range: 760.000, 1200.000, 5
  0.296456256810E+03
 -0.149761277860E+01
  0.317871039240E-02
 -0.318476867010E-05
  0.157208190040E-08
 -0.306913690560E-12
""","""\
name: reference function on ITS-90
type: K
temperature units: °C
emf units: mV
range: -270.000, 0.000, 10
  0.000000000000E+00
  0.394501280250E-01
  0.236223735980E-04
 -0.328589067840E-06
 -0.499048287770E-08
 -0.675090591730E-10
 -0.574103274280E-12
 -0.310888728940E-14
 -0.104516093650E-16
 -0.198892668780E-19
 -0.163226974860E-22
range: 0.000, 1372.000, 9
 -0.176004136860E-01
  0.389212049750E-01
  0.185587700320E-04
 -0.994575928740E-07
  0.318409457190E-09
 -0.560728448890E-12
  0.560750590590E-15
 -0.320207200030E-18
  0.971511471520E-22
 -0.121047212750E-25
exponential:
 a0 =  0.118597600000E+00
 a1 = -0.118343200000E-03
 a2 =  0.126968600000E+03
""","""\
name: reference function on ITS-90
type: N
temperature units: °C
emf units: mV
range: -270.000, 0.000, 8
  0.000000000000E+00
  0.261591059620E-01
  0.109574842280E-04
 -0.938411115540E-07
 -0.464120397590E-10
 -0.263033577160E-11
 -0.226534380030E-13
 -0.760893007910E-16
 -0.934196678350E-19
range: 0., 1300., 10
  0.000000000000E+00
  0.259293946010E-01
  0.157101418800E-04
  0.438256272370E-07
 -0.252611697940E-09
  0.643118193390E-12
 -0.100634715190E-14
  0.997453389920E-18
 -0.608632456070E-21
  0.208492293390E-24
 -0.306821961510E-28
""","""\
name: reference function on ITS-90
type: R
temperature units: °C
emf units: mV
range: -50.000, 1064.180, 9
  0.000000000000E+00
  0.528961729765E-02
  0.139166589782E-04
 -0.238855693017E-07
  0.356916001063E-10
 -0.462347666298E-13
  0.500777441034E-16
 -0.373105886191E-19
  0.157716482367E-22
 -0.281038625251E-26
range: 1064.180, 1664.500, 5
  0.295157925316E+01
 -0.252061251332E-02
  0.159564501865E-04
 -0.764085947576E-08
  0.205305291024E-11
 -0.293359668173E-15
range: 1664.5, 1768.1, 4
  0.152232118209E+03
 -0.268819888545E+00
  0.171280280471E-03
 -0.345895706453E-07
 -0.934633971046E-14
""","""\
name: reference function on ITS-90
type: S
temperature units: °C
emf units: mV
range: -50.000, 1064.180, 8
  0.000000000000E+00
  0.540313308631E-02
  0.125934289740E-04
 -0.232477968689E-07
  0.322028823036E-10
 -0.331465196389E-13
  0.255744251786E-16
 -0.125068871393E-19
  0.271443176145E-23
range: 1064.180, 1664.500, 4
  0.132900444085E+01
  0.334509311344E-02
  0.654805192818E-05
 -0.164856259209E-08
  0.129989605174E-13
range: 1664.5, 1768.1, 4
  0.146628232636E+03
 -0.258430516752E+00
  0.163693574641E-03
 -0.330439046987E-07
 -0.943223690612E-14
""","""\
name: reference function on ITS-90
type: T
temperature units: °C
emf units: mV
range: -270.000, 0.000, 14
  0.000000000000E+00
  0.387481063640E-01
  0.441944343470E-04
  0.118443231050E-06
  0.200329735540E-07
  0.901380195590E-09
  0.226511565930E-10
  0.360711542050E-12
  0.384939398830E-14
  0.282135219250E-16
  0.142515947790E-18
  0.487686622860E-21
  0.107955392700E-23
  0.139450270620E-26
  0.797951539270E-30
range: 0.000, 400.000, 8
  0.000000000000E+00
  0.387481063640E-01
  0.332922278800E-04
  0.206182434040E-06
 -0.218822568460E-08
  0.109968809280E-10
 -0.308157587720E-13
  0.454791352900E-16
 -0.275129016730E-19
"""]

output = '''\
"""
This module contains thermocouple reference functions for types
B,E,J,K,N,R,S,T.

You can access the lookup table objects like so:
    typeK = <this module>.thermocouples['K']

This module contains the NIST ITS-90 thermocouple reference functions.

Disclaimers
-----------
(Author) I make no warranties as to the accuracy of this module, and shall
        not be liable for any damage that may result from errors or omissions.
(NIST) The National Institute of Standards and Technology (NIST) uses its
        best efforts to produce a Database of high quality and to verify that
        the data contained therein have been selected on the basis of sound
        scientific judgement. However, NIST makes no warranties to that effect,
        and NIST shall not be liable for any damage that may result from errors
        or omissions in the Database.

(Note: This module is generated code from create_tables_NIST.py.)
"""

__copyright__ = "public domain"

import numpy as np
from .function_types import Thermocouple_Reference, Polynomial_Gaussian_Piecewise_Function

thermocouples = {

'''

for text in nist_texts:
    ttype, maker = Thermocouple_from_NIST_text(text)
    output += "'"+ttype+"': "+maker+",\n\n"

output += '''}

#end of module'''

print(output)



""" Some strings for describing temperature and voltage units. """
__copyright__ = "public domain"

import numpy as np
import sys

# declare names of the valid temperature and voltage units
Tunits_short = {
    'K': 'K',
    'C': '\xb0C', # python3: use unicode degree symbol
    'R': '\xb0R',
    'F': '\xb0F',
  } if sys.version_info[0] >= 3 else {
    'K': 'K',
    'C': 'degC', # python2: use ascii
    'R': 'degR',
    'F': 'degF',
  }
Tunits_long = {
    'K': 'kelvin',
    'C': 'degrees Celsius',
    'R': 'degrees Rankine',
    'F': 'degrees Fahrenheit',
  }

Vunits_short = {
    'V':  'V',
    'mV': 'mV',
    'uV': '\u03BCV', # python3: use unicode mu symbol
  } if sys.version_info[0] >= 3 else {
    'V':  'V',
    'mV': 'mV',
    'uV': 'uV', # python2: use ascii
  }
Vunits_long = {
    'V':  'volts',
    'mV': 'millivolts',
    'uV': 'microvolts',
  }

# declare unit conversion matrices for temperatures
Tunits_mat = {
'C': {'from': {
            'K': np.array([[1.,-273.15],[0.,1.]]),
            'C': np.array([[1.,0.],[0.,1.]]),
            'R': np.array([[5./9.,-273.15],[0.,1.]]),
            'F': np.array([[5./9.,-32.*5./9.],[0.,1.]]),
              },
      'to':   {
            'K': np.array([[1.,273.15],[0.,1.]]),
            'C': np.array([[1.,0.],[0.,1.]]),
            'R': np.array([[9./5.,491.67],[0.,1.]]),
            'F': np.array([[9./5.,32.],[0.,1.]]),
              },
    },
  }

#Tunits_mat_to_K = {
    #'K': np.array([[1.,0.],[0.,1.]]),
    #'C': np.array([[1.,273.15],[0.,1.]]),
    #'R': np.array([[5./9.,0.],[0.,1.]]),
    #'F': np.array([[5./9.,459.67*5./9.],[0.,1.]]),
  #}
#Tunits_mat_from_K = {
    #'K': np.array([[1.,0.],[0.,1.]]),
    #'C': np.array([[1.,-273.15],[0.,1.]]),
    #'R': np.array([[9./5.,0.],[0.,1.]]),
    #'F': np.array([[9./5.,-459.67],[0.,1.]]),
  #}


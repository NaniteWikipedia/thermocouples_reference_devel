"""
Python module for computing thermocouple emf values from temperatures.
This module just contains the generic thermocouple clas and helper
functions.
"""

__author__    = "User:Nanite @ wikipedia"
__copyright__ = "public domain"

import numpy as np
from .units import *

# scipy.optimize will be imported when needed.
optimize = None
def ensure_import_optimize():
    global optimize
    if optimize == None:
        try:
            import scipy.optimize as optimize
        except ImportError:
            raise ImportError("Inverse lookup requires scipy.optimize module. Please install SciPy.")

class Polynomial_Gaussian_Piecewise_Function(object):
    """\
    Piecewise mathematical function of polynomials plus gaussian, used for
    thermocouple reference.
    
    Main methods:
     func(T)          # compute the function
     func.__call__(T) # synonym for func(T)
     func.inverse(F)  # perform inverse lookup
    
    The raw function parameters are stored in .table. The structure of .table
    is a list of tuples giving the different segments of the piecewise function,
    formatted as:
        (minimum T, maximum T, polynomial coefs array, exponential coefs list)
    The polynomial coefs array is in the order of np.polyval(), i.e., starting
    with the highest power and ending with zeroth power (offset).
    Exponential coefs are used as  ec[0] * np.exp(ec[1] * (T - ec[2])**2), or
    may be None in which case only the polynomial is used.
    
    The appropriate temperature and voltage units to use when calling these
    functions are stored in .Tunits and .Vunits, respectively.
    
    .source and .calibration are strings containing information about where the
    function data comes from, and how it is calibrated.
    """
    def __init__(self, table, Tunits, Vunits, source="", calibration=""):
        self.table       = table
        self.Tunits      = Tunits
        self.Vunits      = Vunits
        self.source      = source
        self.calibration = calibration
        
        # check table
        lastmax = table[0][0]
        for tmin,tmax,pc,ec in table:
            if not tmin <= tmax:
                raise ValueError("Temperature limits must be in ascending order.")
            if tmin != lastmax:
                raise ValueError("Pieces' limits must be contiguous.")
            lastmax = tmax
    
    @property
    def minT(self):
        return self.table[0][0]
    @property
    def maxT(self):
        return self.table[-1][1]
    
    def __repr__(self):
        return "<piecewise polynomial+gaussian, domain {} to {} in {}, output in {}; {} calibrated, from {}>".format(
            self.minT, self.maxT,
            Tunits_short[self.Tunits], Vunits_short[self.Vunits],
            self.calibration, self.source)
    
    def __call__(self,T,derivative=0,out_of_range="raise"):
        """\
        Calculate reference function at given temperature.

        Parameters
        ----------
        T : array_like
            Temperature or array of temperatures.
        derivative: integer
            Use this parameter to evaluate the functional derivative of the emf
            function at a given temperature. Default is derivative=0 (no derivative).
        out_of_range: string, optional
            Determines behaviour for out of range temperatures.
            "raise": raises an ValueError exception. (default)
            "nan":   values replaced by nans.
            "extrapolate": extrapolates from closest range. Do not trust this!
        
        Returns
        -------
        emf : array_like
            computed emf function
        """
        
        if out_of_range not in ["raise", "nan", "extrapolate"]:
            raise ValueError("invalid out_of_range parameter",out_of_range)

        T = np.array(T, copy=False, order='A')
        emf_choices = [None]

        # We go through the table, determining the selector which is used
        # to choose which piece of the piecewise function to use.
        # selector = 0 where T is underrange,
        # selector = 1 where T is in first range,
        #  ...
        # selector = N where T is in last (Nth) range,
        # selector = N+1 where T is overrange.
        tmin = self.minT
        selector = (T >= tmin)*1
        for tmin, tmax, coefs, ec in self.table:
            selector += (T > tmax)
            # Here we go ahead and compute emf values using all ranges.
            #   this is simple but perhaps a bit inefficient.
            emf = np.polyval(np.polyder(coefs, derivative),T)
            
            if ec:
                # Type K thermocouple has this annoying exponential addition term,
                # corresponding to a little bump at 127 Celsius.
                dT = T - ec[2]
                gauss = ec[0] * np.exp(ec[1] * dT**2)
                if derivative == 0:
                    emf += gauss
                elif derivative == 1:
                    emf += 2. * ec[1] * gauss * dT
                elif derivative == 2:
                    emf += 2. * ec[1] * gauss * (2. * ec[1] * dT**2 + 1.)
                elif derivative == 3:
                    emf += 4. * ec[1] * ec[1] * gauss * dT * (2. * ec[1] * dT**2 + 3.)
                else:
                    raise ValueError("sorry, derivatives > 3 not supported for this type.")
            emf_choices.append(emf)
        emf_choices.append(None)
        
        if out_of_range == "nan":
            emf_choices[0] = T*np.nan
            emf_choices[-1] = emf_choices[0]
        else:
            emf_choices[0]  = emf_choices[1]
            emf_choices[-1] = emf_choices[-2]

        if out_of_range == "raise":
            unders = selector <= 0
            overs = selector > len(self.table)
            if np.any(unders) or np.any(overs):
                u_temps = np.extract(unders,T)
                o_temps = np.extract(overs,T)
                if u_temps.size == 0: u_temps = None
                if o_temps.size == 0: o_temps = None
                msg = "Temperatures ("+Tunits_short[self.Tunits]+") under or over range:"
                raise ValueError(msg, u_temps, o_temps)

        return np.choose(selector, emf_choices)

    def inverse(self,V,Tstart=None,Vtol=1e-6):
        """
        Find the temperature corresponding to a given voltage, via zero-finding.
        
        Parameters
        ----------
        V: float
            Measured voltage (in appropriate units) goes here.
        Tstart: float
            Suggested starting temperature for search. If not provided, uses midpoint
            of range. You can speed the search convergence by providing a
            good starting guess here. This function anyway tries its best to converge.
        Vtol: float
            Desired absolute tolerance of voltage value.
        
        Returns
        -------
        T: float
            Temperature T, such that func(T) = V
            Note that the result is checked before returning: if the solution
            would have |func(T) - V| > Vtol, an exception is raised instead.

        Note on implementation
        ----------------------
        First this method tries to use scipy.optimize.newton;
        failing that, it uses scipy.optimize.brentq.
        
        For non-monotonic emf functions this may fail. Among the standard thermocouple
        functions such only occurs with type B in the range 0-50 degC. Anyway, since
        these thermocouples are generally used for higher temperature you should be
        able to avoid that situation.
        
        This function requires scipy to be installed. Upon the first instance
        of calling this function, it attempts to import scipy.optimize.
        This function might take a few milliseconds per call. If that's too
        slow for you, you may want to roll your own solution.
        """

        ensure_import_optimize()
        
        if Tstart == None:
            Tstart = 0.5*(self.minT + self.maxT)
        V = float(V)
        
        # We use "extrapolate" here to allow Newton search to play outside the
        # allowed range in the hope that it returns later on.
        fun0 = lambda T: self(T,out_of_range="extrapolate") - V
        fun1 = lambda T: self(T,derivative=1,out_of_range="extrapolate")
        fun2 = lambda T: self(T,derivative=2,out_of_range="extrapolate")
        try:
            # Try newton's method first.
            T = optimize.newton(fun0, Tstart, fprime=fun1, fprime2=fun2, tol=Vtol)
            if abs(self(T)-V) > Vtol:
                raise ValueError
        except:
            # Any problems (range error, convergence, whatever), then try brentq
            #
            # FIXME: we are assuming the emf function is monotonic here.
            # for type B thermocouple that's not quite true, and so with that type
            # we could raise exception if someone asks for a low V value.
            # (anyway they should not ask such a thing, as type B is meant
            #  for high temperatures.)
            try:
                T = optimize.brentq(fun0, self.minT, self.maxT)
            except ValueError as e:
                if e.args == ("f(a) and f(b) must have different signs",):
                    raise ValueError("Voltage not within in allowed range.")
                else:
                    raise
            if not abs(self(T,out_of_range="nan") - V) <= Vtol:
                raise ValueError("Did not converge within tolerance.")
        
        return T



def doc_emf(uT, uV):
    Tlong = Tunits_long[uT]
    Tshort = Tunits_short[uT]
    Vlong = Vunits_long[uV]
    Vshort = Vunits_short[uV]
    """ Decorator for generating emf_XXX docstrings """
    def fn(f):
        Tref_default = f.__defaults__[0]
        f.__doc__ = """\
        Compute electromotive force for given thermocouple measurement junction
        temperature and given reference junctions temperature.
        
        This method uses {} temperature units and {}.

        Parameters
        ----------
        T : array_like
            Temperature or array of temperatures (in {}).
        Tref : float, optional
            Reference junctions' temperature (in {}),
            defaults to {}.
            If derivative != 0, Tref is irrelevant.
        derivative : integer, optional
            Use this parameter to evaluate the functional derivative of
            the emf function at a given temperature.
            defaults to derivative=0 (no derivative).
        out_of_range : {{'raise', 'nan', 'extrapolate'}}, optional
            Determines behaviour for out of range temperatures: raise an
            exception, return NaNs, or extrapolate using the nearest
            polynomial. Note - do not trust the extrapolation!
        
        Returns
        -------
        emf : array_like
            computed emfs (in {})
            or, if derivative != 0,
            emf derivative (in {} / {}**derivative)
        """.format(Tlong, Vlong, Tshort, Tshort, Tref_default,
                   Vshort, Vshort, Tshort)
        return f
    return fn

def doc_inverse(uT,uV):
    Tlong = Tunits_long[uT]
    Tshort = Tunits_short[uT]
    Vlong = Vunits_long[uV]
    Vshort = Vunits_short[uV]
    """ Decorator for generating inverse_XXX docstrings """
    def fn(f):
        Tref_default = f.__defaults__[0]
        Vtol_default = f.__defaults__[2]
        f.__doc__ = """\
        Inverse lookup: compute measurement junction temperature for a given
        measured voltage and given reference junctions temperature.
        
        You must have SciPy installed to use this method.
        (see documentation of .func.inverse for more notes on implementation)
        
        This method uses {} temperature units and {}.
        
        Parameters
        ----------
        emf : float
            The measured voltage (in {}).
        Tref : float, optional
            The reference junctions' temperature (in {}).
            This allows you to perform cold-junction compensation. Note that
            Tref = {}, the default, corresponds to the reference junctions
            being at the freezing point of water.
        Tstart : float, optional
            Suggested starting temperature (in {}).
            You can hasten the search convergence by providing a good starting
            guess here. If not provided, the midpoint of the entire temperature
            range will be used.
        Vtol : float, optional
            Tolerance of voltage in search (in {}),
            defaults to {}.
        
        Returns
        -------
        T : float
            Junction temperature (in {}), such that:
              emf == func(T) - func(Tref)    (to within Vtol)
        """.format(Tlong, Vlong, Vshort, Tshort, Tref_default, Tshort,
                   Vshort, Vtol_default, Tshort)
        return f
    return fn


class Thermocouple_Reference(object):
    """
    Thermocouple reference helper object. This object provides practical
    methods for converting between temperatures and measured voltages:
    
    * ``.emf_mVC(T)`` converts to millivolts from known temperature in
      degrees Celsius.
    * ``.inverse_CmV(mV)`` converts to degrees Celsius from known voltage
      in millivolts.
    
    Other temperature units are also supported; the complete table being:
    
    ==================   ==========   ==============
     Temperature unit    EMF lookup   Inverse lookup
    ==================   ==========   ==============
    degrees Celsius      .emf_mVC     .inverse_CmV
    degrees Fahrenheit   .emf_mVF     .inverse_FmV
    kelvins              .emf_mVK     .inverse_KmV
    degrees Rankine      .emf_mVR     .inverse_RmV
    ==================   ==========   ==============
    
    In each case it is possible (and desirable) to pass in the reference
    junction temperature by the keyword argument Tref. In practice, that
    junction is often not at the default water-ice point value, and these
    methods take care of the the cold junction compensation in the correct
    way.

    The attribute .func gives access to the raw lookup function. Note that
    the units of this function are typically (but not always) Celsius, so
    you should check them (via .func.Tunits and .func.Vunits). The object
    in .func also contains information about where its data comes from, and
    how its temperature scale is exactly defined (by conventions of 1948,
    1968, or 1990, etc.); such information is useful for traceability and
    precision applications.
    
    The attribute .type is a short string that reads "Type X" for a
    standard letter-designated thermocouple of type X, or is a brief
    codename for thermocouples which have not been assigned a letter.
    
    The attribute .composition is a longer string containing a nominal
    chemical composition of this thermocouple type, in the format
    "PPPP - NNNN", where PPPP is the nominal positive leg (higher
    Seebeck coefficient) and NNNN is the nominal negative leg (lower
    Seebeck coefficient). However note that, strictly speaking, the
    thermocouple types are *not* defined by their composition, but
    rather by their function. For example, a "type K" thermocouple is
    defined as a thermocouple that follows NIST's type K characteristic
    function to the desired degree. Alloy manufacturers use various
    deliberate impurities (dopants) to tune their materials to match
    the standard curves, and in some cases the wire is made of
    an entirely different material (as is common with cheap "extension-
    grade" wire that only matches the curve over a restricted range).
    """
    
    def __init__(self,func,ttype='',composition=''):
        """
        func is the object that contains the actual function information, and has
        methods __call__, inverse, and attributes .minT, .maxT, .Tunits, .Vunits .
        """
        if func.Tunits != 'C' or func.Vunits != 'mV':
            raise ValueError("Only mV <- deg C functions supported right now.")
        
        self.func        = func
        self._mats_Tunits_to   = Tunits_mat['C']['to']
        self._mats_Tunits_from = Tunits_mat['C']['from']
        self.type        = ttype
        self.composition = composition
    
    def __repr__(self):
        rng = "{:.1f} to {:.1f}".format(self.func.minT,self.func.maxT)
        return "<{} thermocouple reference ({} {})>".format(
                self.type, rng, Tunits_short[self.func.Tunits])
    
    @property
    def minT_C(self):
        imul, iadd = self._mats_Tunits_to  ['C'][0]
        return self.func.minT*imul + iadd
    @property
    def maxT_C(self):
        imul, iadd = self._mats_Tunits_to  ['C'][0]
        return self.func.maxT*imul + iadd

    @doc_emf('C','mV')
    def emf_mVC(self,T,Tref=0.,derivative=0,out_of_range="raise"):
        #mul,add = self._mats_Tunits_from['C'][0]
        #T = T*mul + add
        f_T   = self.func(T   ,derivative=derivative,out_of_range=out_of_range)
        if derivative != 0:
            return f_T # * (mul**derivative)
        #Tref = Tref*mul + add
        f_ref = self.func(Tref,derivative=derivative,out_of_range=out_of_range)
        return f_T - f_ref
    
    @doc_emf('F','mV')
    def emf_mVF(self,T,Tref=32.,derivative=0,out_of_range="raise"):
        mul,add = self._mats_Tunits_from['F'][0]
        T = T*mul + add
        f_T   = self.func(T   ,derivative=derivative,out_of_range=out_of_range)
        if derivative != 0:
            return f_T * (mul**derivative)
        Tref = Tref*mul + add
        f_ref = self.func(Tref,derivative=derivative,out_of_range=out_of_range)
        return f_T - f_ref
    
    @doc_emf('K','mV')
    def emf_mVK(self,T,Tref=273.15,derivative=0,out_of_range="raise"):
        mul, add = self._mats_Tunits_from['K'][0]
        T = T*mul + add
        f_T   = self.func(T   ,derivative=derivative,out_of_range=out_of_range)
        if derivative != 0:
            return f_T * (mul**derivative)
        Tref = Tref*mul + add
        f_ref = self.func(Tref,derivative=derivative,out_of_range=out_of_range)
        return f_T - f_ref
    
    @doc_emf('R','mV')
    def emf_mVR(self,T,Tref=491.67,derivative=0,out_of_range="raise"):
        mul, add = self._mats_Tunits_from['R'][0]
        T = T*mul + add
        f_T   = self.func(T   ,derivative=derivative,out_of_range=out_of_range)
        if derivative != 0:
            return f_T * (mul**derivative)
        Tref = Tref*mul + add
        f_ref = self.func(Tref,derivative=derivative,out_of_range=out_of_range)
        return f_T - f_ref
    
    
    @doc_inverse('C','mV')
    def inverse_CmV(self,emf,Tref=0.,Tstart=None,Vtol=1e-6):
        #mul, add = self._mats_Tunits_from['C'][0]
        #Tref = Tref*mul + add
        f_ref = self.func(Tref)
        if Tstart != None: Tstart = Tstart*mul + add
        T = self.func.inverse(emf+f_ref,
                    Tstart=Tstart, Vtol=Vtol)
        #imul, iadd = self._mats_Tunits_to['C'][0]
        #T = T*imul + iadd
        return T
    
    @doc_inverse('F','mV')
    def inverse_FmV(self,emf,Tref=32.,Tstart=None,Vtol=1e-6):
        mul, add = self._mats_Tunits_from['F'][0]
        Tref = Tref*mul + add
        f_ref = self.func(Tref)
        if Tstart != None: Tstart = Tstart*mul + add
        T = self.func.inverse(emf+f_ref,
                    Tstart=Tstart, Vtol=Vtol)
        imul, iadd = self._mats_Tunits_to['F'][0]
        T = T*imul + iadd
        return T
    
    @doc_inverse('K','mV')
    def inverse_KmV(self,emf,Tref=273.15,Tstart=None,Vtol=1e-6):
        mul, add = self._mats_Tunits_from['K'][0]
        Tref = Tref*mul + add
        f_ref = self.func(Tref)
        if Tstart != None: Tstart = Tstart*mul + add
        T = self.func.inverse(emf+f_ref,
                    Tstart=Tstart, Vtol=Vtol)
        imul, iadd = self._mats_Tunits_to['K'][0]
        T = T*imul + iadd
        return T
    
    @doc_inverse('R','mV')
    def inverse_RmV(self,emf,Tref=491.67,Tstart=None,Vtol=1e-6):
        mul, add = self._mats_Tunits_from['R'][0]
        Tref = Tref*mul + add
        f_ref = self.func(Tref)
        if Tstart != None: Tstart = Tstart*mul + add
        T = self.func.inverse(emf+f_ref,
                    Tstart=Tstart, Vtol=Vtol)
        imul, iadd = self._mats_Tunits_to['R'][0]
        T = T*imul + iadd
        return T

#end of module
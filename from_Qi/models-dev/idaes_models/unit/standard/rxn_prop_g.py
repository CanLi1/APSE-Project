# -*- coding: utf-8 -*-
"""
Demonstration property package for the Gibbs rector unit model

Contains the basic set of properties needed to run Gibbs reactor model

Model assumes only a single gas phase with ideal behaviour for calcuating
Gibbs free energy and heat capacities/enthalpies.

@author: alee
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Andrew Lee"
__version__ = "0.1"

#import Pyomo Stuff
from pyomo.environ import *

# List of all chemical species included in package
comp = ['H2','N2','O2','CH4','CO','CO2','H2O','NH3']

# List of all chemical elements that constitute the chemical species 
elem = ['H', 'N', 'O', 'C']

# Elemental composition of all species
elem_comp = {'H2'  : {'H':2, 'N':0, 'O':0, 'C':0},
             'N2'  : {'H':0, 'N':2, 'O':0, 'C':0},
             'O2'  : {'H':0, 'N':0, 'O':2, 'C':0},
             'CH4' : {'H':4, 'N':0, 'O':0, 'C':1},
             'CO'  : {'H':0, 'N':0, 'O':1, 'C':1},
             'CO2' : {'H':0, 'N':0, 'O':2, 'C':1},
             'H2O' : {'H':2, 'N':0, 'O':1, 'C':0},
             'NH3' : {'H':3, 'N':1, 'O':0, 'C':0}
}

# DIPPR ideal gas heat capacity correlation parameters
cpigdip = {'H2'  : {1: 27.617, 2: 9.5600, 3: 2466.0, 4: 3.76000, 5: 567.60},
           'N2'  : {1: 29.105, 2: 8.6149, 3: 1701.6, 4: 0.10347, 5: 909.79},
           'O2'  : {1: 29.103, 2: 10.040, 3: 2526.5, 4: 9.35600, 5: 1153.8},
           'CH4' : {1: 33.298, 2: 79.933, 3: 2086.9, 4: 41.6020, 5: 991.96},
           'CO'  : {1: 29.108, 2: 8.7730, 3: 3085.1, 4: 8.45530, 5: 1538.2},
           'CO2' : {1: 29.370, 2: 34.540, 3: 1428.0, 4: 26.4000, 5: 588.00},
           'H2O' : {1: 33.363, 2: 26.790, 3: 2610.5, 4: 8.89600, 5: 1169.0},
           'NH3' : {1: 33.427, 2: 48.980, 3: 2036.0, 4: 22.5600, 5: 882.00}
}
class PropPack():

    """
    This package contains the necessary property calcualtions to
    demonstrate the basic Gibbs reactor model.
    
    Properties supported:
        - elemental compositions
        - pure component molar Gibbs free energy
    
    """

    def __call__(self, blk, indx=0):
        '''
            Callable method for Block construction
        '''

        self.pyomo_block = blk

        self._make_params(blk)
        self._make_vars_main(blk)

        self._make_constraints_main(blk)

        return blk

    def _make_params(self, blk):
        # Create package parameters
        ''' This section is for parameters needed for the property models.'''
        
        # Component list - a list of component identifiers
        blk.comp = comp
        
        # Component elemental composition
        blk.elem_comp = elem_comp  

        # Gas constant
        blk.R = 8.314   # J/mol.K

        # Mixture heat capacity
        blk.cp = 29.07 # J/mol.K

        # Thermodynamic reference state            
        blk.Pref= Param(within=PositiveReals, default=101325,
                      doc='Reference pressure [Pa]')
        blk.Tref= Param(within=PositiveReals, default=303.15,
                      doc='Reference temperature [K]')
        
    def _make_vars_main(self, blk):
        # Create state variables
        ''' This section contains the state variables to be used to calucate
            the properties, and needs to be linked to the state properties in
            the unit model.'''
            
        blk.P = Var(domain=Reals, initialize=101325,
                    doc='State pressure [Pa]')
        blk.T = Var(domain=Reals, initialize=303.15,
                    doc='State temperature [K]')
        blk.y = Var(blk.comp, domain=Reals, initialize=0,
                    doc='State component mole fractions')

        # ---------------------------------------------------------------------
        # Create standard property variables
        ''' This section is for the creating property variables to be
            calcuated that will be seen the unit model. Generally, these 
            will come from the list of standard properties, and should use
            standard names where possible.'''
        # Pure component vapour heat capacities
        blk.cp_vap_pc = Var(blk.comp, domain=Reals, initialize=1.0,
                        doc="Pure component vapour heat capacities [J/mol.K]")
                        
        # Pure component vapour enthalpies
        blk.h_vap_pc = Var(blk.comp, domain=Reals, initialize=1.0,
                        doc="Pure component vapour enthalpies [J/mol]")

        # Mixture molar enthalpy
        blk.h_mix = Var(domain=Reals, initialize=0,
                        doc='Mixture specific entahlpy [J/mol]')
        
        # Pure component Gibbs free energy
        blk.g_pc = Var(blk.comp, domain=Reals, initialize=0.0, 
                      doc="Pure component Gibbs free energy J/mol")

    def _make_constraints_main(self, blk):
        # Create standard constraints
        ''' This section creates the necessary constraints for calculating
            the standard property values.
            All calcuations assume ideal gas behaviour
        '''
        # Pure component heat capacities
        def rule_cp_vap_pc(b, j):
            return b.cp_vap_pc[j] == cpigdip[j][1] \
                                    + cpigdip[j][2]*((cpigdip[j][3]/b.T) \
                                        / sinh(cpigdip[j][3]/b.T))**2 \
                                    + cpigdip[j][4]*((cpigdip[j][5]/b.T) \
                                        / cosh(cpigdip[j][5]/b.T))**2
        blk.eq_cp_vap_pc = Constraint(blk.comp, rule=rule_cp_vap_pc)
        
        # Pure component molar enthalpy
        def rule_h_pc(b, j):
            return b.h_vap_pc[j] == (cpigdip[j][1]*b.T \
                        + cpigdip[j][2]*cpigdip[j][3]/tanh(cpigdip[j][3]/b.T) \
                        - cpigdip[j][4]*cpigdip[j][5] \
                            * tanh(cpigdip[j][5]/b.T)) \
                        - (cpigdip[j][1]*b.Tref + cpigdip[j][2] \
                            * cpigdip[j][3]/tanh(cpigdip[j][3]/b.Tref) \
                        - cpigdip[j][4]*cpigdip[j][5]*tanh(cpigdip[j][5] \
                            / b.Tref))
        blk.eq_h_pc = Constraint(blk.comp, rule=rule_h_pc)
        
        # Mixture specific enthalpy
        blk.eq_h_mix = Constraint(expr = blk.h_mix \
                                            == sum(blk.y[j]*blk.h_vap_pc[j] \
                                                    for j in blk.comp))

        # Pure component Gibbs free energy
        # Assume constant cp for simplicity and vapour phase only
        def rule_g_pc(b, j):
            return b.g_pc[j] == b.cp_vap_pc[j]*(b.T - b.Tref) \
                                - b.T*b.cp_vap_pc[j]*log(b.T/b.Tref) \
                                + b.T*b.R*log(b.P/b.Pref)
        blk.eq_g_pc = Constraint(blk.comp, rule=rule_g_pc)

def _initialize(blk, T=None, P=None, y=None):
    ''' Initialisation routine for property package'''
        
    if P == None:
        blk.P = 101325
    else:
        blk.P = P
        
    if T == None:
        blk.T = 298.15
    else:
        blk.T = T

    if y == None:
        for j in blk.comp:
            blk.y[j] = 1/len(blk.comp)
    else:
        for j in blk.comp:
            blk.y[j] = y[j]

    for j in blk.comp:
        blk.cp_vap_pc[j] = cpigdip[j][1] \
                            + cpigdip[j][2]*((cpigdip[j][3]/value(blk.T)) \
                                / sinh(cpigdip[j][3]/value(blk.T)))**2 \
                            + cpigdip[j][4]*((cpigdip[j][5]/value(blk.T)) \
                                / cosh(cpigdip[j][5]/value(blk.T)))**2     
        
        blk.h_vap_pc[j] = (cpigdip[j][1]*value(blk.T) \
                        + cpigdip[j][2]*cpigdip[j][3]/tanh(cpigdip[j][3] \
                            / value(blk.T)) \
                        - cpigdip[j][4]*cpigdip[j][5] \
                            * tanh(cpigdip[j][5]/value(blk.T))) \
                        - (cpigdip[j][1]*value(blk.Tref) + cpigdip[j][2] \
                            * cpigdip[j][3]/tanh(cpigdip[j][3] \
                            / value(blk.Tref)) \
                        - cpigdip[j][4]*cpigdip[j][5]*tanh(cpigdip[j][5] \
                            / value(blk.Tref)))
        
        blk.g_pc[j] = value(blk.cp_vap_pc[j])*(value(blk.T) - value(blk.Tref))\
                        - value(blk.T)*value(blk.cp_vap_pc[j]) \
                            * log(value(blk.T)/value(blk.Tref)) \
                        + value(blk.T)*value(blk.R)*log(value(blk.P) \
                            / value(blk.Pref))

    blk.h_mix = sum(value(blk.y[j])*value(blk.h_vap_pc[j]) for j in blk.comp)

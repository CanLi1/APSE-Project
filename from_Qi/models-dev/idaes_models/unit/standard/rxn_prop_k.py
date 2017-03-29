# -*- coding: utf-8 -*-
"""
Demonstration property package for unit rector models

Contains the basic set of properties needed to run reactor models except 
Gibbs energy minimisation unit model.

@author: alee
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Andrew Lee"
__version__ = "0.1"

#import Pyomo Stuff
from pyomo.environ import *

comp = ['a','b','c','d','e','f']
rxn_idx = [1,2,3]  

class PropPack():

    """
    Example property package for reactions

    This package contains the necessary property calcualtions to
    demonstrte the basic unit reactor models.

    System involeves six components (a,b,c,d,e and f) involved in
    three reactions (labled 1, 2 and 3). Reactions equations are:

    a + 2b <-> c + d
    a + 2c <-> 2e
    a + b  <-> f

    Reactions are assumed to be aqueous and only a liquid phase is
    considered.
    
    Properties supported:
        - stoichiometric coefficients
        - rate of reaction
            - rate coefficients (forward and reverse)
            - equilibrium coefficients
        - heats of reaction
        - specific enthalpy of the fluid mixture
    
    """

    def __call__(self, blk, indx=0):
        '''
            Callable method for Block construction
        '''

        self.pyomo_block = blk

        self._make_params(blk)
        self._make_vars_main(blk)
        self._make_vars_sec(blk)

        self._make_constraints_main(blk)
        self._make_constraints_sec(blk)

        return blk

    def _make_params(self, blk):
        # Create package parameters
        ''' This section is for parameters needed for the property models.'''
        
        # Component list - a list of component identifiers
        blk.comp = comp
        
        # Reaction indices - a list of identifiers for each reaction
        blk.rxn_idx = rxn_idx  
        
        # Stoichiometirc coefficients
        '''Stoichiometric coefficient for each component in each reaction'''
        blk.stoic = {
            (1,'a'):-1, (1,'b'):-2, (1,'c'):1, (1,'d'):1, (1,'e'):0, (1,'f'):0,
            (2,'a'):-1, (2,'b'):0, (2,'c'):-2, (2,'d'):0, (2,'e'):2, (2,'f'):0,
            (3,'a'):-1, (3,'b'):-1, (3,'c'):0, (3,'d'):0, (3,'e'):0, (3,'f'):1
        }

        # Gas constant
        blk.R = 8.314   # J/mol.K

        # Mixture heat capacity
        blk.cp = 75.327 # J/mol.K

        # Thermodynamic reference state            
        blk.Pref= Param(within=PositiveReals, mutable=True, default=101325,
                      doc='Reference pressure [Pa]')
        blk.Tref= Param(within=PositiveReals, mutable=True, default=303.15,
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
            
        # Reaction rate
        blk.rate_rxn = Var(blk.rxn_idx, domain=Reals, initialize=0,
                         doc='Normalised Rate of Reaction [mol/m^3.s]')

        blk.dH_rxn = Var(blk.rxn_idx, domain=Reals, initialize=0,
                         doc='Heats of Reaction [J/mol]')

        blk.h_mix = Var(domain=Reals, initialize=0,
                        doc='Mixture specific entahlpy [J/mol]')

    def _make_vars_sec(self, blk):
        # Create package specific variables
        ''' This section is for creating variables that are specific to the
            property package and are for mainly internal use. These do not
            need to follow any standards are they are not expected to be used
            by the units models in most cases.'''
            
        # Reaction reate coefficients
        blk.kf = Var(blk.rxn_idx, doc='Rate coefficient for forward reaction')
        blk.kr = Var(blk.rxn_idx, doc='Rate coefficient for reverse reaction')
        blk.K_eq = Var(blk.rxn_idx, initialize=1,
                                   doc='Equilibrium coefficient')

    def _make_constraints_main(self, blk):
        # Create standard constraints
        ''' This section creates the necessary constraints for calculating
            the standard property values.'''
 
        # Heats of reaction
        ''' For this example, heats of reaction are constants'''
        def rule_dH_rxn(b, i):
            if i == 1:
                return b.dH_rxn[i] == 60000
            elif i == 2:
                return b.dH_rxn[i] == 50000
            else:
                return b.dH_rxn[i] == 80000
        blk.eq_dH_rxn = Constraint(blk.rxn_idx, rule=rule_dH_rxn)
        
        # Reaction rates
        ''' Standard rate expressions for each reaction, with a forward and
            reverse reaction term. 
            Reaction rates are defined as rate of generation, and normalised
            such that the reaction rate for a given component is equal to the
            normalised rate multiplied by the component stoichiometric
            coefficient.'''
        def rule_rate_expr(b, j):
            if j == 0:
                return Constraint.Skip
            else:
                if j == 1:
                    return b.rate_rxn[j] == b.kf[j]*(b.y['a'])*(b.y['b']**2) \
                                            - b.kr[j]*b.y['c']*b.y['d']
                elif j == 2:
                    return b.rate_rxn[j] == b.kf[j]*(b.y['a'])*(b.y['c']**2) \
                                            - b.kr[j]*b.y['e']
                else:
                    return b.rate_rxn[j] == b.kf[j]*(b.y['a'])*(b.y['b']) \
                                            - b.kr[j]*b.y['f']
        blk.eq_rate_rxn = Constraint(blk.rxn_idx, rule=rule_rate_expr)
        
        # Mixture specific enthalpy
        ''' The mixture enthalpy is assumed to be equal to that of pure
            water in the liquid state, with a constant heat capacity'''
        blk.eq_h_mix = Constraint(expr = blk.h_mix == blk.cp* \
                                                        (blk.T - blk.Tref))

    def _make_constraints_sec(self, blk):
        # Create package specific constraints
        ''' This section is for any additional supporting constraints
            required by the property package. These will mainly be
            for calcuating the package specfic variable values.'''

        # Forward rate constants
        ''' Arhenius expression for rate coefficients'''
        def rule_kf(b, i):
            '''if i == 1:
                return b.kf[i] == 17700*exp(-12000/(b.R*b.T))
            elif i == 2:
                return b.kf[i] == 1490*exp(-7000/(b.R*b.T))
            else:
                return b.kf[i] == 26500*exp(-13000/(b.R*b.T))'''
            if i == 1:
                return b.kf[i] == 141.54
            elif i == 2:
                return b.kf[i] == 88.46
            else:
                return b.kf[i] == 139.45
        blk.eq_kf = Constraint(blk.rxn_idx, rule=rule_kf)

        # Reverse rate constants
        ''' Reverse reaction rates coefficients in terms of forward
            coefficients and equilibrium coefficient'''
        def rule_kr(b, i):
            return b.kr[i]*b.K_eq[i] == b.kf[i]
        blk.eq_kr = Constraint(blk.rxn_idx, rule=rule_kr)

        # Equilibrium coefficients
        ''' Equilibrium coefficients as a function of temperature using the
            van't Hoff equation'''
        def rule_K_eq(b, i):
            '''if i == 1:
                return log(b.K_eq[i]) - log(20) == -(b.dH_rxn[i]/b.R)* \
                                                    (1/b.T - 1/298.15)
            elif i == 2:
                return log(b.K_eq[i]) - log(5) == -(b.dH_rxn[i]/b.R)* \
                                                    (1/b.T - 1/298.15)
            else:
                return log(b.K_eq[i]) - log(10) == -(b.dH_rxn[i]/b.R)* \
                                                    (1/b.T - 1/298.15)'''
            if i == 1:
                return b.K_eq[i] == 20
            elif i == 2:
                return b.K_eq[i] == 5
            else:
                return b.K_eq[i] == 10
        blk.eq_K_eq = Constraint(blk.rxn_idx, rule=rule_K_eq)

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

    for i in blk.rxn_idx:
        if i == 1:
            blk.dH_rxn[i] = 60000
            blk.kf[i] = 141.54
            blk.K_eq[i] = 20
        elif i == 2:
            blk.dH_rxn[i] = 50000
            blk.kf[i] = 88.46
            blk.K_eq[i] = 5
        else:
            blk.dH_rxn[i] = 80000
            blk.kf[i] = 139.45
            blk.K_eq[i] = 10
        blk.kr[i] = value(blk.kf[i])/value(blk.K_eq[i])

    for j in blk.rxn_idx:
        if j == 1:
            return blk.rate_rxn[j] == value(blk.kf[j])*(value(blk.y['a']))\
                                            * (value(blk.y['b'])**2) \
                                        - value(blk.kr[j])* value(blk.y['c']) \
                                            * value(blk.y['d'])
        elif j == 2:
            return blk.rate_rxn[j] == value(blk.kf[j])*(value(blk.y['a']))\
                                            * (value(blk.y['c'])**2) \
                                        - value(blk.kr[j])*value(blk.y['e'])
        else:
            return blk.rate_rxn[j] == value(blk.kf[j])*(value(blk.y['a']))\
                                            * (value(blk.y['b'])) \
                                        - value(blk.kr[j])*value(blk.y['f'])
            
    blk.h_mix = value(blk.cp)*(value(blk.T) - value(blk.Tref))
        
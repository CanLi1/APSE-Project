# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 14:24:35 2017

@author: alee
"""
# Chages the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Andrew Lee"
__version__ = "0.1"

# Import Pyomo
from pyomo.environ import *
from pyomo.opt import SolverStatus

# Import IDAES cores
import idaes_models.core.unit_model as unit

comp = ['H2O']

class PropPack(unit.UnitModel):
    """
    This package provedes the necessary constraints for a general cubic
    equation of state.
    """
    def build(self, *args, **kwargs):
        ''' Callable method for Block construction.
        '''
        blk = self

        # Build model
        self._make_params()
        self._make_vars()
        self._make_constraints()
        
        return blk

    def __init__(self, *args, **kwargs):
        """
        Create a property package object.
        """
        # Call base class constructor
        unit.UnitModel.__init__(self, *args, **kwargs)

    def _make_params(self):
        ''' This section contains a number of parameters or constants
            required by the model.
        '''
        blk = self

        # Gas constant
        blk.R = Param(default=8.314, doc='Gas constant [J/mol.K]')
        blk.Tref = Param(default=298.15, doc='Reference Temperature [K]')

    def _make_vars(self):
        """
        Create property variables.
        """
        blk = self
        
        # Set component list
        blk.comp = comp

        # State variables
        blk.T = Var(initialize=298.15)
        blk.P = Var(initialize=101325)
        blk.y = Var(blk.comp, initialize=1)

        # Gas Phase Properties
        blk.h_mix = Var(domain = Reals)
        blk.rho_mix = Var(domain = Reals)

    def _make_constraints(self):
        """
        Create property constraints.
        """
        blk = self

        def rule_eq_t1(b):
            return b.h_mix == 75.2652*(blk.T-blk.Tref)
        blk.eq_t1 = Constraint(rule=rule_eq_t1)

        blk.rho_mix.fix(985.9393497324021)

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

comp = ['CO2','H2O','N2']

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
        # Extract build arguments
        self.prop_list = kwargs.pop("prop_list",{"V_vap",
                                                 "D_vap",
                                                 "k_vap",
                                                 "rho_vap",
                                                 "mu_vap",
                                                 "MW_vap",
                                                 "cp_vap",
                                                 "h_vap"})

        # Call base class constructor
        unit.UnitModel.__init__(self, *args, **kwargs)

    def _make_params(self):
        ''' This section contains a number of parameters or constants
            required by the model.
        '''
        blk = self

        # Gas constant
        blk.R = Param(default=8.314, doc='Gas constant [J/mol.K]')
        
        blk.Tref = Param(default=298.15,
                    doc='Thermodynamic Reference Temperature [K]')

    def _make_vars(self):
        """
        Create property variables.
        """
        blk = self
        
        # Set component list
        blk.comp = comp
        
        blk.T = Var(initialize=298.15)
        blk.P = Var(initialize=101325)
        blk.y = Var(blk.comp, initialize=1)

        # Gas Phase Properties
        if "V_vap" in blk.prop_list:
            blk.V_vap = Var(domain=Reals)
        if "D_vap" in blk.prop_list:
            blk.D_vap = Var(blk.comp, domain=Reals)
        if "k_vap" in blk.prop_list:
            blk.k_vap = Var(domain = Reals)
        if "rho_vap" in blk.prop_list:
            blk.rho_vap = Var(domain=Reals)
        if "mu_vap" in blk.prop_list:
            blk.mu_vap = Var(domain=Reals)
        if "MW_vap" in blk.prop_list or "rho_vap" in blk.prop_list:
            blk.MW_vap = Var(domain=Reals)
        if "cp_vap" in blk.prop_list or "h_vap" in blk.prop_list:
            blk.cp_vap = Var(domain = Reals)
        if "h_vap" in blk.prop_list:
            blk.h_vap = Var(domain = Reals)

    def _make_constraints(self):
        """
        Create property constraints.
        """
        blk = self

        if "V_vap" in blk.prop_list:
            blk.eq_q1 = Constraint(expr = blk.P*blk.V_vap == blk.R*blk.T)
        if "MW_vap" in blk.prop_list or "rho_vap" in blk.prop_list:
            blk.eq_q2 = Constraint(expr = blk.MW_vap == (blk.y['CO2']*44 \
                                + blk.y['H2O']*18 + blk.y['N2']*28)*1e-3)
        if "rho_vap" in blk.prop_list:
            blk.eq_q3 = Constraint(expr = blk.rho_vap == blk.MW_vap*blk.P \
                                                        / (blk.R*blk.T))
        if "D_vap" in blk.prop_list:
            def rule_eq_q4(b, j):
                if j == 'CO2':
                    return b.D_vap[j] == 1e-4*((0.1593-0.1282*(b.P*1e-5-1.4) \
                            + 0.001*(b.T-333.15) \
                            + 0.0964*((b.P*1e-5-1.4)**2) \
                            - 0.0006921*((b.P*1e-5-1.4)*(b.T-333.15)) \
                            - 3.3532e-06*(b.T-333.15)**2)*b.y['H2O'] \
                            / (b.y['H2O']+b.y['N2']) \
                            + (0.1495-0.1204*(b.P*1e-5-1.4) \
                            + 0.0008896*(b.T-333.15) \
                            + 0.0906*((b.P*1e-5-1.4)**2) \
                            - 0.0005857*(b.P*1e-5-1.4)*(b.T-333.15) \
                            - 3.559e-06*(b.T-333.15)**2)*b.y['N2'] \
                            / (b.y['H2O']+b.y['N2']))
                elif j == 'H2O':
                    return b.D_vap[j] == 1e-4*((0.1593-0.1282*(b.P*1e-5-1.4) \
                            + 0.001*(b.T-333.15) \
                            + 0.0964*((b.P*1e-5-1.4)**2) \
                            - 0.0006921*((b.P*1e-5-1.4)*(b.T-333.15)) \
                            - 3.3532e-06*(b.T-333.15)**2)*b.y['CO2'] \
                            / (b.y['CO2']+b.y['N2']) \
                            + (0.2165-0.1743*(b.P*1e-5-1.4) \
                            + 0.001377*(b.T-333.15) \
                            + 0.13109*((b.P*1e-5-1.4)**2) \
                            - 0.0009115*(b.P*1e-5-1.4)*(b.T-333.15) \
                            - 4.8394e-06*(b.T-333.15)**2)*b.y['N2'] \
                            / (b.y['CO2']+b.y['N2']))
                else:
                    return b.D_vap[j] == 1e-4*((0.1495-0.1204*(b.P*1e-5-1.4) \
                            + 0.0008896*(b.T-333.15) \
                            + 0.0906*((b.P*1e-5-1.4)**2) \
                            - 0.0005857*(b.P*1e-5-1.4)*(b.T-333.15) \
                            - 3.559e-06*(b.T-333.15)**2)*b.y['CO2'] \
                            / (b.y['H2O']+b.y['CO2']) \
                            + (0.2165-0.1743*(b.P*1e-5-1.4) \
                            + 0.001377*(b.T-333.15) \
                            + 0.13109*((b.P*1e-5-1.4)**2) \
                            - 0.0009115*(b.P*1e-5-1.4)*(b.T-333.15) \
                            - 4.8394e-06*(b.T-333.15)**2)*b.y['H2O'] \
                            / (b.y['H2O']+b.y['CO2']))
            blk.eq_q4 = Constraint(blk.comp, rule=rule_eq_q4)
        if "k_vap" in blk.prop_list:
            blk.k_vap.fix(0.0280596)
        if "mu_vap" in blk.prop_list:
            blk.mu_vap.fix(1.92403e-5)
        if "cp_vap" in blk.prop_list or "h_vap" in blk.prop_list:
            blk.cp_vap.fix(30.0912)
        if "h_vap" in blk.prop_list:
            blk.eq_q5 = Constraint(expr = blk.h_vap \
                                        == blk.cp_vap*(blk.T-blk.Tref))

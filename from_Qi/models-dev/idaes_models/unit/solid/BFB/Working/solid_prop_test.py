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

comp = ['H2O','Car']
gcomp = ['CO2','H2O','N2']
rxn_idx = [1,2]

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
        self.prop_list = kwargs.pop("prop_list",{"ap",
                                                 "cp_sol",
                                                 "dp",
                                                 "e_mf",
                                                 "k_sol",
                                                 "rho_sol",
                                                 "v_mf",
                                                 "F",
                                                 "r",
                                                 "h_sol"})

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
        blk.scomp = comp
        blk.gcomp = gcomp
        blk.rxn_idx = rxn_idx

        # State variables
        blk.T = Var(initialize=298.15)
        blk.P = Var(initialize=101325)
        blk.y = Var(blk.gcomp, initialize=1)
        blk.w = Var(blk.scomp, initialize=1)

        # Solid Phase Properties
        if "ap" in blk.prop_list:
            blk.ap = Var()
        if "cp_sol" in blk.prop_list or "h_sol" in blk.prop_list:
            blk.cp_sol = Var()
        if "dp" in blk.prop_list or "ap" in blk.prop_list:
            blk.dp = Var()
        if "e_mf" in blk.prop_list:
            blk.e_mf = Var()
        if "k_sol" in blk.prop_list:
            blk.k_sol = Var()
        if "rho_sol" in blk.prop_list:
            blk.rho_sol = Var()
        if "v_mf" in blk.prop_list:
            blk.v_mf = Var()
        if "F" in blk.prop_list:
            blk.F = Var(domain=Reals)
        if "r" in blk.prop_list or "h_sol" in blk.prop_list:
            blk.dH = Var(blk.rxn_idx)
        if "r" in blk.prop_list:
            blk.A = Var(blk.rxn_idx)
            blk.dS = Var(blk.rxn_idx)
            blk.E = Var(blk.rxn_idx)
            blk.nv = Var()
            blk.k = Var(blk.rxn_idx, domain=Reals, initialize=1)
            blk.K = Var(blk.rxn_idx, domain=Reals, initialize=1)
            blk.r = Var(blk.rxn_idx, domain=Reals, initialize=0)
            blk.rg = Var(blk.gcomp, domain=Reals, initialize=0)
            blk.rs = Var(blk.scomp, domain=Reals, initialize=0)
        if "h_sol" in blk.prop_list:
            blk.h_sol = Var(within=Reals, initialize=1)
            blk.cpgcs = Var(blk.gcomp, domain = Reals)

    def _make_constraints(self):
        """
        Create property constraints.
        """
        blk = self

        if "ap" in blk.prop_list:
            blk.ap = 90.49773755656109
            blk.eq_ap = Constraint(expr = blk.ap*blk.rho_sol*blk.dp == 6)
        if "cp_sol" in blk.prop_list or "h_sol" in blk.prop_list:
            blk.cp_sol.fix(1130)
        if "dp" in blk.prop_list or "ap" in blk.prop_list:
            blk.dp.fix(1.5e-4)
        if "e_mf" in blk.prop_list:
            blk.e_mf.fix(0.5)
        if "k_sol" in blk.prop_list:
            blk.k_sol.fix(1.36)
        if "rho_sol" in blk.prop_list:
            blk.rho_sol.fix(442)
        if "v_mf" in blk.prop_list:
            blk.v_mf.fix(0.0085)
        if "F" in blk.prop_list:
            blk.F.fix(0)
        if "r" in blk.prop_list or "h_sol" in blk.prop_list:
            blk.dH[1].fix(-72580.3)
            blk.dH[2].fix(-109691)
        if "r" in blk.prop_list:
            blk.A[1].fix(0.17583)
            blk.A[2].fix(141.944)
            blk.dS[1].fix(-141.425)
            blk.dS[2].fix(-281.255)
            blk.E[1].fix(29622.8)
            blk.E[2].fix(27522.5)
            blk.nv.fix(1900.46)

            def rule_eq_r1(b, i):
                return b.k[i] == b.A[i]*b.T*exp(-b.E[i]/(b.R*b.T))
            blk.eq_r1 = Constraint(blk.rxn_idx, rule=rule_eq_r1)

            def rule_eq_r2(b, i):
                return b.K[i]*b.P == exp(-b.dH[i]/(b.R*b.T) + b.dS[i]/b.R)
            blk.eq_r2 = Constraint(blk.rxn_idx, rule=rule_eq_r2)

            def rule_eq_r3(b, i):
                if i == 1:
                    return b.r[i] == b.k[i]*((b.P*b.y['H2O']) \
                                - (b.w['H2O']*b.rho_sol/b.K[i]))
                else:
                    return b.r[i] == b.k[i]*(((1-2*(b.w['Car']*b.rho_sol/b.nv))**2) \
                            * (b.P*b.y['CO2']) \
                            - ((b.w['Car']*b.rho_sol/b.nv) \
                            * ((b.w['Car']*b.rho_sol/b.nv))/b.K[i]))
            blk.eq_r3 = Constraint(blk.rxn_idx, rule=rule_eq_r3)

            def rule_eq_r4(b, j):
                if j == 'CO2':
                    return b.rg[j] == b.nv*b.r[2]
                elif j == 'H2O':
                    return b.rg[j] == b.r[1]
                else:
                    return b.rg[j] == 0
            blk.eq_r4 = Constraint(blk.gcomp, rule=rule_eq_r4)

            def rule_eq_r5(b, j):
                if j == 'H2O':
                    return b.rs[j] == b.r[1]
                else:
                    return b.rs[j] == b.nv*b.r[2]
            blk.eq_r5 = Constraint(blk.scomp, rule=rule_eq_r5)
            
        if "h_sol" in blk.prop_list:
            blk.cpgcs['CO2'].fix(38.7)
            blk.cpgcs['H2O'].fix(34.1)
            blk.cpgcs['N2'].fix(29.2)

            blk.eq_s2 = Constraint(expr = blk.h_sol == ((blk.w['H2O']) \
                        * (blk.cpgcs['H2O']*(blk.T-blk.Tref) + blk.dH[1]) \
                        + blk.w['Car']*(blk.cpgcs['CO2']*(blk.T-blk.Tref) \
                            + blk.dH[2])) \
                        + blk.cp_sol*(blk.T-blk.Tref))

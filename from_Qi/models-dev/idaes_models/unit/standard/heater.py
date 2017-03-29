"""
Standard IDAES heater models.
"""
from __future__ import division
from __future__ import print_function

__author__ = "John Eslick, Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.core.base.connector import IndexedConnector

import idaes_models.core.unit_model as unit

class Heater(unit.UnitModel):
    """
    7. Heater/Cooler Unit Class
    
    Model has 1 degree of freedom (with fixed inlet stream).
    Usually set: deltaT, Out_T, or Q.
    """
    
    _required_properties = [
        ("h_vap", "vapor_phase enthalpy")]
    
    def __init__(self, comp, prop_dict):
        """
        Create a heater/cooler object with the component set comp.
        
        Args
        ====
        comp: a list of components (should be in fixed order)
        prop_lib: a default external function library for properties
        prop_lib: dictionary with unit model property name as key and
            a list of [library file name, alternate name or none]
        """
        unit.UnitModel.__init__(self, prop_dict=prop_dict)
        self.comp = comp #component set
        
    def _make_vars(self, b):
        """
        Create the heater models block variables
        """
        # Inlet stream
        b.In_F = Var(domain=NonNegativeReals, initialize=1.0,
                     doc="Inlet total molar flow rate (mol/s)")
        b.In_P = Var(domain=NonNegativeReals, initialize=101325,
                     doc="Inlet pressure (Pa)")
        b.In_T = Var(domain=NonNegativeReals, initialize=273.15,
                     doc="Inlet temperature (K)")
        b.In_y = Var(self.comp, domain=NonNegativeReals, initialize=1.0,
                     doc="Inlet component mole fractions")
        #Outlet stream
        b.Out_F = Var(domain=NonNegativeReals, initialize=1.0, 
                      doc="Outlet molar flowrate (mol/s)")
        b.Out_P = Var(domain=NonNegativeReals, initialize=101325,
                      doc="Outlet Pressure (Pa)")
        b.Out_T = Var(domain=NonNegativeReals, initialize=273.15,
                      doc="Outlet temperature (K)")
        b.Out_y = Var(self.comp, domain=NonNegativeReals, initialize=1,
                      doc="Outlet component mole fractions")
        # Performance Variables
        b.deltaT = Var(domain=Reals, initialize=0,
                       doc="Temperature Difference (K)")
        b.Q = Var(domain=Reals, initialize=0, doc="Heat duty (J/s)")
                
    def _make_energy_balance(self, b):
        """
        Make energy balance equations.
        """
        # Links to Property Package
        p_h_vap = self.prop_dict["h_vap"]
        
        # Argument Lists for Enthalpy Calls (for neatness)
        arg_in = [b.In_T, b.In_P]
        arg_out = [b.Out_T, b.Out_P]
        for j in self.comp:
            arg_in.append(b.In_y[j])
            arg_out.append(b.Out_y[j])
            
        #Energy balance equation
        b.eq_ebal = Constraint(
            expr = b.In_F*p_h_vap(*arg_in) + b.Q == b.Out_F*p_h_vap(*arg_out))
        # DeltaT
        b.eq_deltaT = Constraint(expr=b.deltaT == b.Out_T - b.In_T)
            
    def _make_mass_balance(self, b):
        """
        Make total and component mass balance equations.
        """
        # Total Balance
        b.eq_tbal = Constraint(expr=b.In_F == b.Out_F)
        # Component Balances
        def rule_mbal(blk, j):
            return blk.In_y[j] == blk.Out_y[j]
        b.eq_cbal = Constraint(self.comp, rule=rule_mbal)
        
    def _make_pressure(self, b):
        """
        No pressure drop over unit: Pin = Pout
        """
        b.eq_P = Constraint(expr = b.Out_P == b.In_P)
        
    def _make_connectors(self, b):
        """
        Make inlet and outlet stream connectors.
        """
        b.Inlet = self.one_phase_fluid_connector(
            b.In_F, b.In_T, b.In_P, b.In_y)
        b.Outlet = self.one_phase_fluid_connector(
            b.Out_F, b.Out_T, b.Out_P, b.Out_y)
    
    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block
        Args
        ====
        b: the pyomo block to set up.
        index: the index if creating an indexed block
        """
        self._make_vars(b)
        self._make_energy_balance(b)
        self._make_mass_balance(b)
        self._make_pressure(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

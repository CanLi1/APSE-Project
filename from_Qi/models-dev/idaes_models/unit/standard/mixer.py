"""
Standard IDAES mixer models.
"""

from __future__ import division
from __future__ import print_function

__author__ = "John Eslick, Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.core.base.connector import IndexedConnector

import idaes_models.core.unit_model as unit

class Mixer(unit.UnitModel):
    """
    1. Mixer Unit Class
    
    Model has 0 degrees of freedom with fixed inlet streams.
    """
    _required_properties = [
        ("h_vap", "vapor_phase enthalpy")]
    
    def __init__(self, comp, prop_dict, NI=2):
        """
        Create a mixer object with the component set comp and NI inlets.
        Args
        ====
        comp: a Pyomo Set containing a list of components
        NI: the number of inlets to the mixer (a positive integer, default=2)
        """
        unit.UnitModel.__init__(self, prop_dict=prop_dict) # Call base class constructor
        self.comp = comp # Component set
        self.NI = list(range(1,NI+1)) #Convert to set of inlet connectors
        #Reference to the main Pyomo model for property calls

    def _make_vars(self, b):
        """
        Make mixer variables
        """
        # Inlet steams: self.NI = [1,2, ... n]
        b.In_F = Var(self.NI, domain=NonNegativeReals, initialize=1.0,
                     doc="Inlet total molar flow rates (mol/s)")
        b.In_P = Var(self.NI, domain=NonNegativeReals, initialize=101325,
                     doc="Inlet pressures (Pa)")
        b.In_T = Var(self.NI, domain=NonNegativeReals, initialize=273.15,
                     doc="Inlet temperatures (K)")
        b.In_y = Var(self.NI, self.comp, domain=NonNegativeReals, 
                     initialize=1.0, doc="Inlet component mole fractions")
        # Outlet stream       
        b.Out_F = Var(domain=NonNegativeReals, initialize=2.0, 
                      doc="Outlet molar flowrate (mol/s)")
        b.Out_P = Var(domain=NonNegativeReals, initialize=101325,
                      doc="Outlet Pressure (Pa)")
        b.Out_T = Var(domain=NonNegativeReals, initialize=273.15,
                      doc="Outlet temperature (K)")
        b.Out_y = Var(self.comp, domain=NonNegativeReals, initialize=1,
                      doc="Outlet component mole fractions")
        # Intermediate Variables
        b.minP = Var(self.NI, domain=NonNegativeReals, initialize=101325,
                     doc="Variable to find minimum inlet pressure (Pa)")
    
    def _make_pressure(self, b):
        """
        Smoth min to find minimum pressure and calculate outlet pressure.
        """
        b.eps = Param(within=PositiveReals, default=1e-4, mutable=True,
                      doc="Smooth minimum smoothing factor")
        # Find minimum inlet pressure by
        #  sequential smooth minimum functions for each inlet
        def rule_minP(blk, i):
            if i == 1:
                return blk.minP[i] == blk.In_P[i]
            else:
                return blk.minP[i] == blk.minP[i-1] -\
                    0.5*(blk.minP[i-1]-blk.In_P[i] +\
                    ((blk.minP[i-1]-blk.In_P[i]+blk.eps)**2)**0.5)
        b.eq_minP = Constraint(self.NI, rule=rule_minP)
        b.eq_P = Constraint(expr = b.Out_P == b.minP[len(self.NI)])
        
    def _make_mass_balance(self, b):
        """
        Make total and component mass balances.
        """
        # Total Balance
        b.eq_tbal = Constraint(
            expr=sum(b.In_F[i] for i in self.NI) == b.Out_F)
        # Component Balances
        def rule_mbal(blk, j):
            return sum(blk.In_F[i]*blk.In_y[i,j] for i in self.NI) ==\
                    blk.Out_F*blk.Out_y[j]
        b.eq_cbal = Constraint(self.comp, rule=rule_mbal)
        
    def _make_connectors(self, b):
        """
        Make inlet and outlet connectors.  Inlet connector set is 
        defined by self.NI
        """
        b.Inlet = Connector(self.NI)
        self.one_phase_fluid_connector(
            b.In_F, b.In_T, b.In_P, b.In_y, index=self.NI, c=b.Inlet)
        b.Outlet = self.one_phase_fluid_connector(
            b.Out_F, b.Out_T, b.Out_P, b.Out_y)
            
    def _make_energy_balance(self, b):
        """
        Make the energy balance equations
        """
        # Links to Property Package
        p_h_vap = self.prop_dict["h_vap"]
        # Argument Lists for Property Calls
        arg_in = {}
        for i in self.NI:
            arg_in[i] = [b.In_T[i], b.In_P[i]]
            for j in self.comp:
                arg_in[i].append(b.In_y[i,j])        
        arg_out = [b.Out_T, b.Out_P]
        for j in self.comp:
            arg_out.append(b.Out_y[j])
        # Energy balance
        b.eq_ebal = Constraint(
            expr=sum(b.In_F[i]*p_h_vap(*arg_in[i]) for i in self.NI) ==
            b.Out_F*p_h_vap(*arg_out))
        
    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block
        Args
        ====
        b: the pyomo block to set up.
        index: the index if creating an indexed block
        """
        self._make_vars(b)
        self._make_mass_balance(b)
        self._make_pressure(b)
        self._make_energy_balance(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

"""
Standard IDAES separator models. (unreviewed)

Can be a splitter if the frac variables for every stream are equal to one
another. i.e., frac[i,j] = x[i] for all i steams for all j components; 
"""
from __future__ import division
from __future__ import print_function

__author__ = "Tony Burgard"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.core.base.connector import IndexedConnector

import idaes_models.core.unit_model as unit

class Separator(unit.UnitModel):
    """
    1. Separator unit class
    """
    def __init__(self, comp, Nout=2, main_model=None):
        """
        Create a separator object with the component set comp and Nout outlets.
        Args
        ====
        comp: a Pyomo Set containing a list of components
        Nout: the number of outlets from the mixer 
            (a positive integer, default=2)
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp
        
        self.Nout = list(range(1,Nout+1)) #Convert to set of outlet connectors
        #Reference to the main Pyomo model for property calls
        self.main_model=main_model
        
    def _make_vars(self, b):
        """
        Make separator variables
        """
        # Inlet stream       
        b.In_F = Var(domain=NonNegativeReals, initialize=2.0, 
                      doc="Inlet molar flowrate (mol/s)")
        b.In_P = Var(domain=NonNegativeReals, initialize=101325,
                      doc="Inlet pressure (Pa)")
        b.In_T = Var(domain=NonNegativeReals, initialize=273.15,
                      doc="Inlet temperature (K)")
        b.In_y = Var(self.comp, domain=NonNegativeReals, initialize=1,
                      doc="Inlet component mole fractions")

        # Outlet steams: self.Nout = [1,2, ... n]
        b.Out_F = Var(self.Nout, domain=NonNegativeReals, initialize=1.0,
                     doc="Outlet total molar flow rates (mol/s)")
        b.Out_P = Var(self.Nout, domain=NonNegativeReals, initialize=101325,
                     doc="Outlet pressures (Pa)")
        b.Out_T = Var(self.Nout, domain=NonNegativeReals, initialize=273.15,
                     doc="Outlet temperatures (K)")
        b.Out_y = Var(self.Nout, self.comp, domain=NonNegativeReals, 
                     initialize=1.0, doc="Outlet component mole fractions")
        b.frac =  Var(self.Nout, self.comp,domain=NonNegativeReals, 
                              initialize=0.5,
                              doc="Mole fraction of each component in "\
                                  "each outlet stream")
        
    def _make_pressure(self, b):
        """
        No pressure drop over unit: Pin = Pout
        """
        b.eq_P = Constraint (self.Nout, noruleinit=True)
        for i in self.Nout:
            b.eq_P.add(i,expr = b.Out_P[i] == b.In_P)
        
    def _make_temperature(self,b):
        """
        Fixing the outlet temperature to be equal to the inlet temperature
        """
        def rule_Temp(blk,i):
            if i > 1:
                return blk.In_T == blk.Out_T[i]
            if i == 1:
                return Constraint.Feasible
        b.eq_fixtemp = Constraint(self.Nout,rule=rule_Temp)    
            
    def _make_mass_balance(self, b):
        """
        Make total and component mass balances.
        """            
        # Total Balance
#        b.eq_tbal = Constraint(
#            expr=sum(b.Out_F[i] for i in self.Nout) == b.In_F)
        
        # Component Balances
        def rule_mbal(blk, j):
            return sum(blk.Out_F[i]*blk.Out_y[i,j] for i in self.Nout) == blk.In_F*blk.In_y[j]
        b.eq_cbal = Constraint(self.comp, rule=rule_mbal)

        def rule_sumy(blk,i):
            return sum(blk.Out_y[i,j] for j in self.comp) == 1
        b.sumy = Constraint(self.Nout, rule=rule_sumy)

        def rule_sumfrac(blk,j):
            return sum(blk.frac[i,j] for i in self.Nout) == 1
        b.sumfrac= Constraint(self.comp, rule=rule_sumfrac)        

        def rule_split(blk,i,j):
            if i > 1:
                return blk.Out_F[i]*blk.Out_y[i,j] == blk.In_F*blk.frac[i,j]*blk.In_y[j]
            if i == 1:
                return Constraint.Feasible
        b.split= Constraint(self.Nout, self.comp, rule=rule_split)         
               
    def _make_connectors(self, b):
        """
        Make inlet and outlet connectors.  Outlet connector set is 
        defined by self.Nout
        """
        b.Inlet = self.one_phase_fluid_connector(b.In_F, b.In_T, b.In_P, b.In_y)
#        b.Outlet = Connector(self.Nout)
#        b.Outlet = self.one_phase_fluid_connector(b.Out_F, b.Out_T, b.Out_P, b.Out_y, index=self.Nout, c=b.Outlet)
            
    def _make_energy_balance(self, b):
        """
        Make the energy balance equations
        """
        # Links to Property Package
        p_h_vap = self.main_model.p_h_vap
        # Argument Lists for Property Calls
        arg_out = {}
        for i in self.Nout:
            arg_out[i] = [b.Out_T[i], b.Out_P[i]]
            for j in self.comp:
                arg_out[i].append(b.Out_y[i,j])        
        arg_in = [b.In_T, b.In_P]
        for j in self.comp:
            arg_in.append(b.In_y[j])
        # Energy balance
        b.eq_ebal = Constraint(
            expr=sum(b.Out_F[i]*p_h_vap(*arg_out[i]) for i in self.Nout) ==
            b.In_F*p_h_vap(*arg_in))
        
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
        self._make_temperature(b)
        self._make_energy_balance(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

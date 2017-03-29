"""
Standard IDAES compressor model.
"""
from __future__ import division
from __future__ import print_function

__author__ = "Jinliang Ma"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.core.base.connector import IndexedConnector

import idaes_models.core.unit_model as unit

class Compressor(unit.UnitModel):
    """
    8. Compressor Unit Class
    
    Model has 2 degree of freedom (with fixed inlet stream).
    Usually set: Ratio_P, Eta (efficiency) and calculate Out_T_Isen, Out_T, Work_Isen, Work
    """
    def __init__(self, comp, main_model=None):
        """
        Create a compressor object with the component set comp.
        
        Args
        ====
        comp: a list of components (should be in fixed order)
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp #component set
        #Reference to the main Pyomo model for property calls
        self.main_model=main_model
        
    def _make_vars(self, b):
        """
        Create model's block variables
        """
        # Inlet stream
        b.In_F = Var(domain=NonNegativeReals, initialize=1.0,
                     doc="Inlet total molar flow rates (mol/s)")
        b.In_P = Var(domain=NonNegativeReals, initialize=101325,
                     doc="Inlet temperatures (Pa)")
        b.In_T = Var(domain=NonNegativeReals, initialize=273.15,
                     doc="Inlet temperatures (K)")
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
        b.Ratio_P = Var(domain=NonNegativeReals, initialize=1,
                       doc="Ratio of outlet to inlet absolute pressure")
        b.Eta = Var(domain=NonNegativeReals, bounds=(0.5,0.99), initialize=0.8, doc="Effciency with respect to isentropic process")
        b.Out_T_Isen = Var(domain=NonNegativeReals, initialize=273.15, doc="Isentropic outlet temperature (K)")
        b.Work = Var(domain=Reals, initialize=0, doc="Work input to compressor")
        b.Work_Isen = Var(domain=Reals, initialize=0, doc="Work input to compressor if isentropic process")
                
    def _make_energy_balance(self, b):
        """
        Make energy balance equations.
        """
        # Links to Property Package
        p_h_vap = self.main_model.p_h_vap
        p_s_vap = self.main_model.p_s_vap
        
        # Argument Lists for Enthalpy and Entropy Calls (for neatness)
        arg_in = [b.In_T, b.In_P]
        arg_out = [b.Out_T, b.Out_P]
        arg_out_isen = [b.Out_T_Isen, b.Out_P]
        for j in self.comp:
            arg_in.append(b.In_y[j])
            arg_out.append(b.Out_y[j])
            arg_out_isen.append(b.Out_y[j])
            
        # Isentropic process, equal entropy
        b.eq_s = Constraint(expr=p_s_vap(*arg_in) == p_s_vap(*arg_out_isen))
        
        # Isentropic work
        b.eq_work_isen = Constraint(expr=b.Work_Isen == b.In_F*(p_h_vap(*arg_out_isen)-p_h_vap(*arg_in)))
        
        # Actual work
        b.eq_work = Constraint(expr=b.Work_Isen == b.Work*b.Eta)
            
        #Energy balance equation
        b.eq_ebal = Constraint(
            expr=b.In_F*p_h_vap(*arg_in) + b.Work == b.Out_F*p_h_vap(*arg_out))
            
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
        Equation to calculate outlet pressure based on pressure ratio
        """
        b.eq_P = Constraint(expr = b.Out_P == b.In_P*b.Ratio_P)
        
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

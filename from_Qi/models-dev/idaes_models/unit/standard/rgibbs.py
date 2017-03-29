"""
Standard IDAES Gibbs reactor model.

Inherits from rstoic and adds Gibbs free energy minimisation,
with supporting properties and elemental balances
"""
from __future__ import division
from __future__ import print_function

__author__ = "Jinliang Ma, Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.unit_model as unit

import rxn_prop_g as prop_mod

class UnitModel(unit.UnitModel):
    """
    Gibbs Reactor Unit Class
    
    This model assume all possible reactions reach equilibrium such that the 
    system partial molar Gibbs free energy is minimized.
    Since some species mole flow rate might be very small, 
    the natural log of the species molar flow rate is used.
    Instead of specifying the system Gibbs free energy as an objective
    function, the equations for zero partial derivatives of the grand function
    with Lagrangian multiple terms with repect to product species mole flow
    rates and the multiples are specified as constraints.
    
    Model has degrees of freedom equal to 2.
    Usually set: - reactor outlet temperature or heat duty,
                 - reactor outlet pressure
    """
    def __init__(self, main_model=None):
        """
        Create a reactor object with the component set comp.
        
        Args
        ====
        comp: a list of components (should be in fixed order)
        """
        unit.UnitModel.__init__(self) #call base class constructor

        #Reference to the main Pyomo model for property calls
        self.main_model=main_model

    def _make_vars(self, b):
        """
        Create the reactor model block variables
        """
        # Create component list set
        b.comp = prop_mod.comp
        b.elem = prop_mod.elem
        
        # Inlet stream
        b.In_F = Var(domain=Reals, initialize=1.0,
                     doc="Inlet total molar flow rate (mol/s)")
        b.In_P = Var(domain=Reals, initialize=101325.0,
                     doc="Inlet pressure (Pa)")
        b.In_T = Var(domain=Reals, initialize=273.15,
                     doc="Inlet temperature (K)")
        b.In_y = Var(b.comp, domain=Reals, initialize=1.0,
                     doc="Inlet component mole fractions")
        # Outlet stream
        b.Out_F = Var(domain=Reals, initialize=1.0, 
                      doc="Outlet molar flowrate (mol/s)")
        b.Out_P = Var(domain=Reals, initialize=101325.0,
                      doc="Outlet pressure (Pa)")
        b.Out_T = Var(domain=Reals, initialize=273.15,
                      doc="Outlet temperature (K)")
        b.Out_y = Var(b.comp, domain=NonNegativeReals, initialize=1.0,
                      doc="Outlet component mole fractions")
        # Performance Variables
        b.Fe = Var(b.elem, domain=Reals, initialize=1.0,
                       doc="Elemental flow rate [mol/s]")
        b.Li = Var(b.elem, domain=Reals, initialize=100,
                      doc="Lagrangian multipliers")
        b.Q = Var(domain=Reals, initialize=0.0, doc="Heat duty (J/s)")

    def _make_prop_block(self, b):
        """ Create Pyomo blocks to contain property calculations
            This overwrites the _make_prop_block method of the base class"""

        # Link to callable class in property package
        prop_rule = prop_mod.PropPack()
        
        # Create property block
        b.prop_in = Block(rule=prop_rule)
        b.prop_out = Block(rule=prop_rule)
        
        # Link state variables
        b.eq_prop_in_T = Constraint(expr = b.In_T == b.prop_in.T)
        b.eq_prop_in_P = Constraint(expr = b.In_P == b.prop_in.P)
        def rule_prop_in_y(b, j):
            return b.In_y[j] == b.prop_in.y[j]
        b.eq_prop_in_y = Constraint(b.comp, rule=rule_prop_in_y)
        
        b.eq_prop_out_T = Constraint(expr = b.Out_T == b.prop_out.T)
        b.eq_prop_out_P = Constraint(expr = b.Out_P == b.prop_out.P)
        def rule_prop_out_y(b, j):
            return b.Out_y[j] == b.prop_out.y[j]
        b.eq_prop_out_y = Constraint(b.comp, rule=rule_prop_out_y)

    def _make_mass_balance(self, b):
        """
        Make total and elemental mass balance equations.
        """
        # Total elemental feed rates
        def rule_elem_feed(blk, e):
            return blk.Fe[e] == blk.In_F*sum(blk.In_y[j] \
                            * blk.prop_in.elem_comp[j][e] for j in blk.comp)
        b.eq_elem_feed = Constraint(b.elem, rule=rule_elem_feed)
        
        # Elemental balances
        def rule_elem_bal(blk, e):
            return blk.Fe[e] == blk.Out_F*sum(blk.Out_y[j] \
                            * blk.prop_out.elem_comp[j][e] for j in blk.comp)
        b.eq_elem_bal = Constraint(b.elem, rule=rule_elem_bal)
        
        # Sum of mole fractions
        b.eq_y_out = Constraint(expr = 1 == sum(b.Out_y[j] for j in b.comp))

    def _make_energy_balance(self, b):
        """
        Make energy balance equations.
        """
        # Energy balance equation
        b.eq_ebal = Constraint(
            expr = 0 == b.In_F*b.prop_in.h_mix - b.Out_F*b.prop_out.h_mix
                            + b.Q)

    def _make_perform(self, b):
        """
        Make performance constraints
        """
        # Use Lagrangian multiple method to derive equations for Out_Fi
        # Use RT*Li as the Lagrangian multiple such that Li is in the
        # similar order of magnitude as log(Yi)
        def rule_lag(blk, j):
            # Use natural log of species mole flow to avoid Pyomo solver
            # warnings of reaching infeasible point
            return 0 == blk.prop_out.g_pc[j] + blk.prop_out.R*blk.Out_T \
                            * (log(blk.Out_y[j]) \
                            + sum(blk.Li[e]*blk.prop_out.elem_comp[j][e] \
                                for e in blk.elem))
        b.eq_lag = Constraint(b.comp, rule=rule_lag)
        
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
        self._make_prop_block(b)
        self._make_mass_balance(b)
        self._make_energy_balance(b)
        self._make_perform(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

    def _initialize(blk, T=None, y=None):
        """ Initialisation procedure for unit model.
        
            Assumes: - inlet conditions are known
                     - reactor pressure is fixed
        """
        
        # Setting initial values for outlet state variables
        if blk.Out_T.fixed == False:        
            if T == None:
                blk.Out_T = value(blk.In_T)
            else:
                blk.T = T

        for j in blk.comp:
            if y == None:
                blk.Out_y[j] = value(blk.In_y[j])
            else:
                blk.Out_y[j] = y[j]

        blk.Out_F = value(blk.In_F)

        # Creating argument lists to pass to property initialisations
        ylist_in = {}
        ylist_out = {}
        for j in blk.comp:
            ylist_in[j] = value(blk.In_y[j])
            ylist_out[j] = value(blk.Out_y[j])

        # Initialize property blocks
        prop_mod._initialize(blk.prop_in,T=value(blk.In_T),
                                         P=value(blk.In_P),y=ylist_in)
        prop_mod._initialize(blk.prop_out,T=value(blk.Out_T),
                                         P=value(blk.Out_P),y=ylist_out)
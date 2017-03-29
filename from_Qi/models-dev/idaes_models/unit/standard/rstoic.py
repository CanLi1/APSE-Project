"""
Standard IDAES stoichiometric reactor model.
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.unit_model as unit

import rxn_prop_k as prop_mod

class UnitModel(unit.UnitModel):
    """
    Stoichiometric Reactor Unit Class
    
    Model has degrees of freedom equal to 2 + number of reaction occuring
    (with fixed inlet stream).
    Usually set: - reactor outlet temperature or heat duty,
                 - reactor outlet pressure
                 - FOR each reaction occuring in reactor:
                     - one extent of reaction for the given reaction OR,
                     - one outlet mole fraction.
    
    Notes:
        - Conversion cannot be defined in a general form, as it is singular
            for species with zero mole fraction in the inlet
        - Yields and selectivities are also difficult to define in a general
            form as they are linked to some key species
        - Users may define these constraints separately if desired
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
        b.rxn_idx = prop_mod.rxn_idx
        
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
        b.x_rxn = Var(b.rxn_idx, domain=Reals, initialize=0.0,
                      doc="Extents of reaction (mol/s)")

        b.Q = Var(domain=Reals, initialize=0.0, doc="Heat duty (J/s)")

    def _make_prop_block(self, b):
        """ Create Pyomo blocks to contain property calculations"""

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
        Make total and component mass balance equations.
        """
        # Component Balances
        def rule_mbal(blk, j):
            return blk.In_F*blk.In_y[j] == blk.Out_F*blk.Out_y[j] \
                                - sum(blk.x_rxn[i]*blk.prop_out.stoic[i,j] for i in blk.prop_out.rxn_idx)
        b.eq_cbal = Constraint(b.comp, rule=rule_mbal)
        
        # Sum of Mole Fractions
        b.eq_sum_y = Constraint(expr = 1 == sum(b.Out_y[j] for j in b.comp))
                
    def _make_energy_balance(self, b):
        """
        Make energy balance equations.
        """
        # Energy balance equation
        b.eq_ebal = Constraint(
            expr = 0 == b.In_F*b.prop_in.h_mix - b.Out_F*b.prop_out.h_mix
                            + sum(b.x_rxn[i]*b.prop_out.dH_rxn[i]
                                        for i in b.prop_out.rxn_idx)
                            + b.Q)
        
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

        y_test = 0
        for j in blk.comp:
            if blk.Out_y[j].fixed == True:
                y_test += 1

        if y_test == 0:
            for j in blk.comp:
                if y == None:
                    blk.Out_y[j] = value(blk.In_y[j])
                else:
                    blk.Out_y[j] = y[j]
        else:
            for j in blk.comp:
                if blk.Out_y[j].fixed == False:
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

        # Initialize remaining unit variables
        if blk.Q.fixed == False:
            blk.Q = -value(blk.In_F)*value(blk.prop_in.h_mix) \
                    + value(blk.Out_F)*value(blk.prop_out.h_mix) \
                    - sum(value(blk.x_rxn[i])*value(blk.prop_out.dH_rxn[i])
                                for i in blk.prop_out.rxn_idx)

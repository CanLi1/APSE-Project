"""
Standard IDAES CSTR reactor model.

Inherits from rstoic and adds rates of reaction, with supporting properties
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.unit_model as unit
import rstoic as rstoic
import rxn_prop_k as prop_mod

class UnitModel(rstoic.UnitModel):
    """
    CSTR Reactor Unit Class
    
    Model has degrees of freedom equal to 3.
    Usually set: - reactor outlet temperature or heat duty,
                 - reactor outlet pressure, 
                 - reactor volume.
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

    def _make_perform(self, b):
        """ Add additional constraints to relate extents of reactions to
            property calculations"""
        
        # Add additional variable for reactor volume
        b.V = Var(domain=Reals, initialize=0.0, doc="Reactor volume (m^3)")
        
        # Add constraint relating extents to rate of reaction
        def rule_x_rxn(blk, i):
            return blk.x_rxn[i] == blk.V*blk.prop_out.rate_rxn[i]
        b.eq_x_rxn = Constraint(b.rxn_idx, rule=rule_x_rxn)

    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block
        Args
        ====
        b: the pyomo block to set up.
        index: the index if creating an indexed block
        """
        # Call the base class
        rstoic.UnitModel.__call__(self, b, index)
        
        # Call additional methods
        self._make_perform(b)
        
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

        # Initialize remaining unit variables
        for i in blk.rxn_idx:
            blk.x_rxn[i] = value(blk.V)*value(blk.prop_out.rate_rxn[i])

        if blk.Q.fixed == False:
            blk.Q = -value(blk.In_F)*value(blk.prop_in.h_mix) \
                    + value(blk.Out_F)*value(blk.prop_out.h_mix) \
                    - sum(value(blk.x_rxn[i])*value(blk.prop_out.dH_rxn[i])
                                for i in blk.prop_out.rxn_idx)
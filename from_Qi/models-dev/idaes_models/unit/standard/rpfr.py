"""
Standard IDAES PFR reactor model.

Base model includes support for extension to include:
    - pressure drop along reactor, deltaP
    - heat duty along reactor length (Q)
    - side inlets along reactor (F_side, y_side and h_side)

By default, all of these variables are fixed to 0, but these can be unfixed as
needed by users and developers of more advanced models.
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.dae import *

import idaes_models.core.unit_model as unit

import rxn_prop_k as prop_mod

class UnitModel(unit.UnitModel):
    """
    PFR Reactor Unit Class
    
    Model has degrees of freedom equal to 2
    Usually set: - reactor length,
                 - reactor area or volume
    """
    def __init__(self, ndp=20, main_model=None):
        """
        Create a reactor object with the component set comp.
        
        Args
        ====
        comp: a list of components (should be in fixed order)
        ndp:  number of discretisation points to use in model
        """
        unit.UnitModel.__init__(self) #call base class constructor

        self.ndp = ndp   # number of points for discretisation
        
        #Reference to the main Pyomo model for property calls
        self.main_model=main_model
        
    def _make_vars(self, b):
        """
        Create the reactor model block variables
        """
        # Create component list set
        b.comp = prop_mod.comp
        b.rxn_idx = prop_mod.rxn_idx

        # Discreitsation domain
        b.l = ContinuousSet(bounds=(0,1))
        
        # Inlet stream
        b.In_F = Var(domain=Reals, initialize=1.0,
                     doc="Inlet total molar flow rate (mol/s)")
        b.In_P = Var(domain=Reals, initialize=101325,
                     doc="Inlet pressure (Pa)")
        b.In_T = Var(domain=Reals, initialize=273.15,
                     doc="Inlet temperature (K)")
        b.In_y = Var(b.comp, domain=Reals, initialize=1.0,
                     doc="Inlet component mole fractions")
        # Outlet stream
        b.Out_F = Var(domain=Reals, initialize=1.0, 
                      doc="Outlet molar flowrate (mol/s)")
        b.Out_P = Var(domain=Reals, initialize=101325,
                      doc="Outlet pressure (Pa)")
        b.Out_T = Var(domain=Reals, initialize=273.15,
                      doc="Outlet temperature (K)")
        b.Out_y = Var(b.comp, domain=Reals, initialize=1,
                      doc="Outlet component mole fractions")
        # Reactor Design and Performance Variables
        b.A = Var(domain=Reals, initialize=1, doc="Reactor area (m^2)")
        b.Lr = Var(domain=Reals, initialize=1, doc="Reactor length (m)")
        b.V = Var(domain=Reals, initialize=1, doc="Reactor volume (m^3)")
        # Discretised Variables
        b.Fy = Var(b.comp, b.l, domain=Reals, initialize=1.0,
                     doc="Component molar flow rate (mol/s)")
        b.Fh = Var(b.l, domain=Reals, initialize=1.0,
                     doc="Enthalpy flow rate (J/s)")
        b.P = Var(b.l, domain=Reals, initialize=101325,
                     doc="Pressure (Pa)")
        b.T = Var(b.l, domain=Reals, initialize=273.15,
                     doc="Temperature (K)")
        b.y = Var(b.comp, b.l, domain=Reals, initialize=1.0,
                     doc="Component mole fractions")
        b.deltaP = Var(b.l, domain=Reals, initialize=0,
                     doc="Pressure drop (Pa/m)")
        b.x_rxn = Var(b.l, b.rxn_idx, domain=Reals, initialize=0,
                      doc="Extents of reaction (mol/m.s)")
        # Expansion feature variables
        b.F_side = Var(b.l, domain=Reals, initialize=0,
                     doc="Side stream total molar flow rate (mol/s)")
        b.y_side = Var(b.comp, b.l, domain=Reals, initialize=1,
                     doc="Side stream component mole fractions")
        b.Fh_side = Var(b.l, domain=Reals, initialize=0,
                     doc="Side stream entahlpy flow rate (J/s)")
        b.Q = Var(b.l, domain=Reals, initialize=0, doc="Heat duty (J/s)")
        # Derivative Variables
        b.dFydl = DerivativeVar(b.Fy, wrt=b.l)
        b.dFhdl = DerivativeVar(b.Fh, wrt=b.l)
        b.dPdl = DerivativeVar(b.P, wrt=b.l)

    def _make_prop_block(self, b):
        """ Create Pyomo blocks to contain property calculations"""

        # Link to callable class in property package
        prop_rule = prop_mod.PropPack()
        
        # Create property block
        b.prop = Block(b.l, rule=prop_rule)
        
        # Link state variables
        def rule_prop_T(b, i):
            return b.T[i] == b.prop[i].T
        b.eq_prop_T = Constraint(b.l, rule=rule_prop_T)
        def rule_prop_P(b, i):
            return b.P[i] == b.prop[i].P
        b.eq_prop_P = Constraint(b.l, rule=rule_prop_P)
        def rule_prop_y(b, i, j):
            return b.y[j,i] == b.prop[i].y[j]
        b.eq_prop_y = Constraint(b.l, b.comp, rule=rule_prop_y)
            
    def _make_mass_balance(self, b):
        """
        Make total and component mass balance equations.
        """
        # Component mass balances
        def rule_mbal(blk, i, j):
            if i == 0:
                return blk.Fy[j,i] == blk.In_F*blk.In_y[j]
            else:
                return 0 == -blk.dFydl[j,i]/blk.Lr \
                                + sum(blk.x_rxn[i,k]*blk.prop[i].stoic[k,j] \
                                                for k in blk.rxn_idx) \
                                + blk.F_side[i]*blk.y_side[j,i]
        b.eq_mbal = Constraint(b.l, b.comp, rule=rule_mbal)
        
        # Mole fractions
        def rule_molefrac(blk, i, j):
            return blk.y[j,i] == blk.Fy[j,i]/sum(blk.Fy[j,i] for j in blk.comp)
        b.eq_molefrac = Constraint(b.l, b.comp, rule=rule_molefrac)

    def _make_energy_balance(self, b):
        """
        Make energy balance equations.
        """
        # Enthalpy flow rate constraint
        def rule_enthflow(blk, i):
            return blk.Fh[i] == sum(blk.Fy[j,i] for j in blk.comp) \
                                    * blk.prop[i].h_mix
        b.eq_enthflow = Constraint(b.l, rule=rule_enthflow)
        
        # Energy balance constraint
        def rule_ebal(blk, i):
            if i == 0:
                return blk.T[i] == blk.In_T
            else:
                return 0 == -blk.dFhdl[i]/blk.Lr + blk.Fh_side[i] + blk.Q[i] \
                        + sum(blk.x_rxn[i,j]*blk.prop[i].dH_rxn[j] \
                                            for j in blk.rxn_idx)
        b.eq_ebal = Constraint(b.l, rule=rule_ebal)

    def _make_rxn_perform(self, b):
        """
        Make reaction performance constraints
        """
        # Extents of reaction
        def rule_x_rxn(blk, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return blk.x_rxn[i,j] == blk.A*blk.prop[i].rate_rxn[j]
        b.eq_x_rxn = Constraint(b.l, b.rxn_idx, rule=rule_x_rxn) 
        
        # Outlet variables
        b.eq_Out_F = Constraint(expr = b.Out_F == \
                                            sum(b.Fy[j,1] for j in b.comp))
        b.eq_Out_T = Constraint(expr = b.Out_T == b.T[1])
        b.eq_Out_P = Constraint(expr = b.Out_P == b.P[1])
        def rule_Out_y(blk, j):
            return blk.Out_y[j] == blk.y[j,1]
        b.eq_Out_y = Constraint(b.comp, rule=rule_Out_y)
        
        # Reactor geometry
        b.eq_geometry = Constraint(expr = b.V == b.Lr*b.A)

    def _make_pressure(self, b):
        """
        Pressure drop constraints
        """
        def rule_deltaP(blk, i):
            if i == 0:
                return blk.P[i] == blk.In_P
            else:
                return blk.dPdl[i] == blk.deltaP[i]*blk.Lr
        b.eq_deltaP = Constraint(b.l, rule=rule_deltaP)

    def _make_transform(self, b):
        b.discretizer = TransformationFactory('dae.finite_difference')
        b.discretizer.apply_to(b,nfe=self.ndp,wrt=b.l,scheme='BACKWARD')
        
        # Fix variables for expansion features
        # This needs to be done after the transformation to fix the internal
        # points.
        b.F_side.fix(0)
        b.y_side.fix(1/len(b.l)) # To avoid potential singularities
        b.Fh_side.fix(0)
        b.Q.fix(0)
        b.deltaP.fix(0)
        
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
        self._make_transform(b)
        self._make_prop_block(b)
        self._make_mass_balance(b)
        self._make_energy_balance(b)
        self._make_pressure(b)
        self._make_rxn_perform(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

    def _initialize(blk, T=None, y=None):
        """ Initialisation procedure for unit model.
        
            Assumes: - inlet conditions are known
        """
        
        # Setting initial values for state variables
        # Assume total flowrate constant
        for i in range(1,len(blk.l)+1):
            if T == None:
                blk.T[blk.l[i]] = value(blk.In_T)
            else:
                blk.T[blk.l[i]] = T

            if i == 1:
                blk.P[blk.l[i]] = value(blk.In_P)
            else:
                blk.P[blk.l[i]] = value(blk.P[blk.l[i-1]]) \
                            - value(blk.deltaP[blk.l[i]]) \
                            * value(blk.Lr)/(len(blk.l)-1)
                
            for j in blk.comp:
                if y == None:
                    blk.y[j,blk.l[i]] = value(blk.In_y[j])
                    blk.Fy[j,blk.l[i]] = value(blk.In_F)*value(blk.In_y[j])
                else:
                    blk.y[j,blk.l[i]] = y[j]
                    blk.Fy[j,blk.l[i]] = value(blk.In_F)*y[j]

        # Initialise outlet variables
        blk.Out_F = value(blk.In_F)
        blk.Out_T = value(blk.T[1])
        blk.Out_P = value(blk.P[1])
        for j in blk.comp:
            blk.Out_y[j] = value(blk.y[j,1])

        # Initialize property blocks
        for i in blk.l:
            # Creating argument lists to pass to property initialisations
            ylist = {}
            for j in blk.comp:
                ylist[j] = value(blk.y[j,i])
            
            prop_mod._initialize(blk.prop[i],T=value(blk.T[i]),
                                         P=value(blk.P[i]),y=ylist)

        # Initialize remaining unit variables

        if blk.A.fixed == False:
            blk.A = value(blk.V)/value(blk.Lr)
        elif blk.V.fixed == False:
            blk.V = value(blk.Lr)*value(blk.A)
        else:
            blk.Lr = value(blk.V)/value(blk.A)

        for i in blk.l:
            blk.Fh[i] = sum(value(blk.Fy[j,i]) for j in blk.comp) \
                            * value(blk.prop[i].h_mix)

            for k in blk.rxn_idx:
                blk.x_rxn[i,k] = value(blk.A)*value(blk.prop[i].rate_rxn[k])

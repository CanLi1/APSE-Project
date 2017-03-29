"""
Standard IDAES flash models. (unreviewed)

"""
from __future__ import division
from __future__ import print_function

__author__ = "Tony Burgard"
__version__ = "1.0.0"

from pyomo.environ import *
from pyomo.core.base.connector import IndexedConnector
#from properties_PR4 import bbEOS,aaEOS,wEOS,uEOS,phiEOS,K,H,Zequations,Zequations2

import idaes_models.core.unit_model as unit

class flash(unit.UnitModel):
    """
    Equilibrium Flash Unit Class
    
    Model has degrees of freedom equal to 2.
    Usually set: - flash outlet temperature or heat duty,
                 - flash outlet pressure
    """
    def __init__(self, comp, Nout=2, main_model=None):
        """
        Create a splitter object with the component set comp and Nout outlets.
        Args
        ====
        comp: a Pyomo Set containing a list of components
        Nout: the number of outlets from the splitter 
            (a positive integer, default=2)
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp
        
        self.Nout = list(range(1,Nout+1)) #Convert to set of inlet connectors
        #Reference to the main Pyomo model for property calls
        self.main_model=main_model
        
    def _make_vars(self, b):
        """
        Make splitter variables
        """
        # Inlet stream       
        b.In_F = Var(domain=NonNegativeReals, initialize=2.0, 
                      doc="Inlet molar flowrate (mol/s)")
        b.In_P = Var(domain=NonNegativeReals, initialize=7*101325,
                      doc="Inlet pressure (Pa)")
        b.In_T = Var(domain=NonNegativeReals, initialize=95,
                      doc="Inlet temperature (K)")

        b.In_y = Var(self.comp, domain=NonNegativeReals, initialize=1/3,
                      doc="Inlet component mole fractions")
        #Outlet stream
        b.Out_F = Var(self.Nout, domain=NonNegativeReals, initialize=1.0, 
                      doc="Outlet molar flowrate (mol/s)")
        b.Out_P = Var(self.Nout, domain=NonNegativeReals, initialize=7*101325,
                      doc="Outlet pressure (Pa)")
        b.Out_T = Var(self.Nout, domain=NonNegativeReals, initialize=105,
                      doc="Outlet temperature (K)")
        b.Out_y = Var(self.Nout, self.comp, domain=NonNegativeReals, 
                      initialize=1/3, doc="Outlet component mole fractions")
        b.Q = Var(domain=Reals, initialize=0, doc="Heat duty (J/s)")
        b.T = Var(domain=NonNegativeReals, initialize=105)
        b.P = Var(domain=NonNegativeReals, initialize=7*101325)

        # Required to relax equilibrium constraint when phase vanishes
        b.sV = Var(domain=NonNegativeReals)
        b.sL = Var(domain=NonNegativeReals)
        b.beta = Var(initialize=1)
        
        b.Alpha = Var(initialize=1)
        b.Beta = Var(initialize=1)
        
        # Will be replaced by C function that calculates Z
        b.ZEOSin = Var(domain=NonNegativeReals,initialize=0.1)
        b.ZEOSout1 = Var(domain=NonNegativeReals,initialize=1)        
        b.ZEOSout2 = Var(domain=NonNegativeReals,initialize=0.1)
                
    def _make_energy_balance(self, b):
        '''
        Make energy balance equations
        '''                                                                   
        b.eq_ebal = Constraint(expr = 0 == b.Q +
            b.In_F*H(b.In_y['O2'],b.In_y['N2'],b.In_y['Ar'],b.In_T,b.In_P,b.ZEOSin) -\
            b.Out_F[1]*H(b.Out_y[1,'O2'],b.Out_y[1,'N2'],b.Out_y[1,'Ar'],b.Out_T[1],b.Out_P[1],b.ZEOSout1) -\
            b.Out_F[2]*H(b.Out_y[2,'O2'],b.Out_y[2,'N2'],b.Out_y[2,'Ar'],b.Out_T[2],b.Out_P[2],b.ZEOSout2))
        
    def _make_mass_balance(self, b):
        """
        Make total and component mass balance equations.
        """
        # Total Balance
        b.eq_tbal = Constraint(
            expr=sum(b.Out_F[i] for i in self.Nout) == b.In_F)
        
        # Will be replaced by a function to calculate Alpha and Beta
        b.Alpha.fix(1.0)
        b.Beta.fix(1.0)
        
        # Component Balances
        def rule_mbal(blk, j):
        #    if j != 'Ar':
                return b.Out_F[1]*blk.Out_y[1,j]*blk.Beta +  b.Out_F[2]*blk.Out_y[2,j]*blk.Alpha== \
                            blk.In_F*blk.In_y[j]
        #    else:
        #        return Constraint.Skip
        b.eq_cbal = Constraint(self.comp, rule=rule_mbal)
        
        def rule_sumy(blk,i):
            return sum(blk.Out_y[i,j] for j in self.comp) == 1
#            return sum(blk.Out_y[i,j]*blk.Out_F[i] for j in self.comp) == blk.Out_F[i]
        b.sumy = Constraint(self.Nout, rule=rule_sumy)

    def _make_flshequations(self, b):
        """
        Make flash equations
        """
        # Outlet streams in temperature equilibrium
        b.EqTRelThrmE = Constraint(expr = b.Out_T[1] == b.Out_T[2])

        # Outlet streams in pressure equilibrium
        b.EqPRelThrmE = Constraint(expr = b.Out_P[1] == b.Out_P[2])

        # Flash temperature sets outlet temperature
        b.EqTRel1Flsh = Constraint(expr = b.Out_T[1] == b.T)

        # Flash pressure sets outlet pressure
        b.EqPRel1Flsh = Constraint(expr = b.Out_P[1] == b.P)
              
        def EqVLEThrm_rule(blk,j):
                return blk.Out_y[1,j] == blk.beta*K(j,blk.Out_y[1,'O2'],blk.Out_y[1,'N2'],blk.Out_y[1,'Ar'],blk.Out_T[1],blk.Out_P[1],blk.ZEOSout1,blk.Out_y[2,'O2'],blk.Out_y[2,'N2'],blk.Out_y[2,'Ar'],blk.Out_T[2],blk.Out_P[2],blk.ZEOSout2)*blk.Out_y[2,j]
        b.EqVLEThrm = Constraint(self.comp,rule=EqVLEThrm_rule)

        # Will be replaced by function to calculate Z
        Zequations(b, self.main_model.ceos_z_liq, self.main_model.ceos_z_vap, self.main_model.ceos_z_liq)
        #Zequations2(b)
        
        b.EqSlackV = Constraint(expr = b.beta - 1 <= b.sV)
        b.EqSlackL = Constraint(expr = -b.sL <= b.beta - 1) 
        
        b.Obj = Objective(expr = b.sL+b.sV,sense=minimize)
        
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
        self._make_flshequations(b)
        self._make_connectors(b)
        self.pyomo_block = b
        return b

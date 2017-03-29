# -*- coding: utf-8 -*-
"""
General Cubic Equation of State model for Pyomo.

This proerty package is designed to be inheritied by a user made class
containing a number of parameters and methods necesary for this model.

The necessary parameters and methods are:
    EoS_ID: an identifier of the desiered Cubic EoS to use
        - Peng-Robinson:        PR
        - Souave-Redlich-Kwong: SRK
    comp: a dictionary of component identifiers (component list)
    Tc: dictionary containing the critical temperatures (K) for all components
    Pc: dictionary containing the critical pressures (Pa) for all components
    omega: dictionary containing the accentricity factors for all components
    kappa: dictionary containing the relevant binary interaction parameters
            for all components

    eq_h_ig_pc: a set of Pyomo constaints indexed by comp returning a value
            for the ideal gas enthalpy of each component [J/mol]
    eq_s_ig_pc: a set of Pyomo constaints indexed by comp returning a value
            for the ideal gas entropy of each component [J/mol.K]
        
@author: alee based on tburgard and lbiegler
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Andrew Lee, Tony Burgard", "Larry Biegler"
__version__ = "0.1"

# Import Pyomo
from pyomo.environ import *
from pyomo.opt import SolverStatus

# Import IDAES cores
import idaes_models.core.unit_model as unit

# -----------------------------------------------------------------------------
# Dictionary of EoS identifiers and parameters
EoS_key = {'PR' : 0,
          'SRK': 1
}
EoS_param = {'PR'  : {'u':2,'w':-1,'omegaA':0.45724},
             'SRK' : {'u':1,'w': 0,'omegaA':0.42748}
}
class PropPack(unit.UnitModel):
    """
    This package provedes the necessary constraints for a general cubic
    equation of state.
    """
    def build(self, *args, **kwargs):
        ''' Callable method for Block construction.
            Args:
            plib = an external library containing external function calls
                    for returning the roots of a cubic EoS.
            phases = flag indicating which phases to include in model.
                    Options V and L ignore potential phase equilibria.
                        - VL: include both liquid and vapour phases
                        - L: consider only liquid phase properties
                        - V: consider only vapour phase properties
            mix_props = flag indicating whether to calculate combined
                    properties for two phase mixtures (requires phases = VL)
                        - 1: calculate mixture properties
                        - 0: do not calculate mixture properties
        '''
        blk = self

        # Set construction flags
        if self.flag_phases != 'VL':
            self.flag_mprop = 0

        # Build model
        self._make_params()
        self._thermo_params()

        self._make_state_vars()
        self._make_external(self.plib)

        self._make_EoS_coeff()

        if self.flag_phases == 'VL':
            self._make_equil()
            self._make_equil_liq()
            self._make_equil_vap()

        self._make_ideal_props()
        self._ideal_pc_props()

        if self.flag_phases != 'V':
            self._make_liq_props()

        if self.flag_phases != 'L':
            self._make_vap_props()

        if self.flag_mprop == 1:
            self._make_thermo_mix()

        self._trans_params()
        if self.flag_phases != 'V':
            self._make_trans_liq()
        if self.flag_phases != 'L':
            self._make_trans_vap()
        if self.flag_mprop == 1:
            self._make_trans_mix()

        return blk

    def __init__(self, *args, **kwargs):
        """
        Create a property package object.
        """
        # Call base class constructor
        unit.UnitModel.__init__(self, *args, **kwargs)

    def _make_params(self):
        ''' This section contains a number of parameters or constants
            required by the model.
        '''
        blk = self

        # Gas constant
        blk.R = Param(default=8.314, doc='Gas constant [J/mol.K]')
        
        # Thermodynamic reference state            
        blk.Pref= Param(within=PositiveReals, default=101325,
                      doc='Reference pressure [Pa]')
        blk.Tref= Param(within=PositiveReals, default=298.15,
                      doc='Reference temperature [K]')
        
        # Smooth IF smoothing parameter
        blk.eps = Param(default=1e-3, mutable=True, doc='Smooth IF parameter')
        
    def _make_state_vars(self):
        # Create state variables
        ''' This section contains the state variables to be used to calucate
            the properties, and needs to be linked to the state properties in
            the unit model.
        '''
        blk = self

        blk.P = Var(domain=Reals, initialize=101325,
                    doc='State pressure [Pa]')
        blk.T = Var(domain=Reals, initialize=303.15,
                    doc='State temperature [K]')

        if blk.flag_phases != 'V':
            blk.x = Var(blk.comp, domain=Reals, initialize=0,
                    doc='Liquid phase component mole fractions')
        if blk.flag_phases != 'L':
            blk.y = Var(blk.comp, domain=Reals, initialize=0,
                    doc='Vapour phase component mole fractions')
        blk.z = Var(blk.comp, domain=Reals, initialize=0,
                    doc='Total mixture component mole fractions')

    def _make_external(self, plib):
        # Set up calls to external procedures
        """ This section sets up the calls to external functions needed to
            perform procedural determination of the compressibility factors.
        """
        blk = self

        # Set up external function calls
        if blk.flag_phases != 'V':
            blk.proc_Z_liq = ExternalFunction(
                library=plib, function="ceos_z_liq")
        if blk.flag_phases != 'L':
            blk.proc_Z_vap = ExternalFunction(
                library=plib, function="ceos_z_vap")

    def _make_EoS_coeff(self):
        # Create coefficients for general cubic EoS
        """ This section uses the user provided parameters to calcuate the
            coefficients of the cubic EoS.
            This section automatically selects the correct form for the
            coefficients based on EoS_ID.
        """
        blk = self

        # Parameters
        # EoS specific parameters
        blk.omegaA = EoS_param[blk.EoS_ID]['omegaA']
        blk.EoS_u = EoS_param[blk.EoS_ID]['u']
        blk.EoS_w = EoS_param[blk.EoS_ID]['w']
        blk.EoS_p = sqrt(blk.EoS_u**2 - 4*blk.EoS_w)

        # Create expressions for coefficients
        def func_fw(b,j):
            if b.EoS_ID == 'PR':
                return 0.37464 + 1.54226*b.omega[j] - 0.26992*b.omega[j]**2
            elif b.EoS_ID == 'SRK':
                return 0.48 + 1.574*b.omega[j] - 0.176*b.omega[j]**2
            else:
                print('Unrecognised EoS_ID')

        def func_b(b,j):
            if b.EoS_ID == 'PR':
                return 0.07780*b.R*b.Tc[j]/b.Pc[j]
            elif b.EoS_ID == 'SRK':
                return 0.08664*b.R*b.Tc[j]/b.Pc[j]
            else:
                print('Unrecognised EoS_ID')
        
        def func_a(b, j):
            return b.omegaA * ((b.R*b.Tc[j])**2/b.Pc[j]) \
                        * ((1+b.fw[j]*(1-sqrt(b.T/b.Tc[j])))**2)

        # Create parameters for constant coefficients
        blk.b = Param(blk.comp, initialize=func_b,
                        doc='Component b coefficient')
        blk.fw = Param(blk.comp, initialize=func_fw,
                        doc='EoS S factor')
        blk.a = Expression(blk.comp, rule=func_a,
                        doc='Component a coefficient')

        # ---------------------------------------------------------------------
        # Variables
        if blk.flag_phases != 'V':
            blk.am_liq = Var(domain=Reals,
                         doc='Liquid phase mixture a coefficient')
            blk.bm_liq = Var(domain=Reals,
                         doc='Liquid phase mixture b coefficient')

        if blk.flag_phases != 'L':
            blk.am_vap = Var(domain=Reals,
                         doc='Vapour phase mixture a coefficient')
            blk.bm_vap = Var(domain=Reals,
                         doc='Vapour phase mixture b coefficient')

        # ---------------------------------------------------------------------
        # Constraints
        # Pressure independent coefficients
        if blk.flag_phases != 'V':
            blk.eq_am_liq = Constraint(expr = blk.am_liq == sum(sum(
                            blk.x[i]*blk.x[j]
                            * sqrt(blk.a[i]*blk.a[j])
                            * (1-blk.kappa[i][j])
                            for j in blk.comp) for i in blk.comp))
            blk.eq_bm_liq = Constraint(expr = blk.bm_liq
                                == sum(blk.x[i]*blk.b[i] for i in blk.comp))

        if blk.flag_phases != 'L':
            blk.eq_am_vap = Constraint(expr = blk.am_vap == sum(sum(
                            blk.y[i]*blk.y[j]
                            * sqrt(blk.a[i]*blk.a[j])
                            * (1-blk.kappa[i][j])
                            for j in blk.comp) for i in blk.comp))
            blk.eq_bm_vap = Constraint(expr = blk.bm_vap
                                == sum(blk.y[i]*blk.b[i] for i in blk.comp))

    def _make_equil(self):
        # Create phase equilibrium variables and constraints
        ''' This section creates the necessary variables and constraints for
            calculating the phase equilibrium in the mixture.
            Fugacity coefficients are defined here as they are needed for the
            equilibrium constraint. Calcuation of the fugacity coefficients
            will be done in the next methods.
        '''
        blk = self

        # Variables
        # Additional pressure variables for smooth behaviour
        blk.Pe = Var(domain=Reals, initialize=101325,
                    doc='Equilibrium pressure [Pa]')

        # Fugacity coefficient
        blk.phi_liq = Var(blk.comp, domain=Reals, initialize=1,
                    doc='Liquid phase fugacity coefficients')
        blk.phi_vap = Var(blk.comp, domain=Reals, initialize=1,
                    doc='Vapour phase fugacity coefficients')

        # Material variables
        blk.L = Var(domain=Reals, initialize=1,
                    doc='Liquid fraction')
        blk.V = Var(domain=Reals, initialize=1,
                    doc='Vapour fraction')

        # Parameters for determining Pe
        blk.beta = Var(domain=Reals, initialize=1,
                        doc='Equilibrium pressure adjustment parameter')
        blk.GL = Var(domain=Reals, initialize=0,
                        doc='Complementarity variable for finding beta')
        blk.GV = Var(domain=Reals, initialize=0,
                        doc='Complementarity variable for finding beta')

        # ---------------------------------------------------------------------
        # Constraints

        # Equilibrium constraint
        def rule_Keq(b, j):
            return b.x[j]*b.phi_liq[j] == b.y[j]*b.phi_vap[j]
        blk.eq_Keq = Constraint(blk.comp, rule=rule_Keq)

        # Material balance equations
        blk.eq_mbal_t = Constraint(expr = 1 == blk.L + blk.V)

        blk.eq_sum_xy = Constraint(expr = 0 == sum(blk.y[j]-blk.x[j]
                                                    for j in blk.comp))

        def rule_mbal_c(b, j):
            return b.z[j] == b.L*b.x[j] + b.V*b.y[j]
        blk.eq_mbal_c = Constraint(blk.comp, rule=rule_mbal_c)

        # Equilibrium state pressure constraint
        blk.eq_Pe = Constraint(expr = blk.Pe*blk.beta == blk.P)

        # Equilibrium pressure adjustment constraint
        blk.eq_beta = Constraint(expr = 0 ==0.5*blk.GV*(blk.beta-1 \
                            + ((blk.beta-1)**2 + (blk.eps)**2)**0.5) \
                            + 0.5*blk.GL*(1-blk.beta \
                            + ((1-blk.beta)**2 + (blk.eps)**2)**0.5))
        blk.eq_GL = Constraint(expr = blk.GL == 0.5*(1-blk.beta + blk.L \
                            - ((blk.L - 1 + blk.beta)**2 + (blk.eps)**2)**0.5))
        blk.eq_GV = Constraint(expr = blk.GV == 0.5*(blk.beta-1 + blk.V \
                            - ((blk.V - blk.beta + 1)**2 + (blk.eps)**2)**0.5))

    def _make_equil_liq(self):
        # Create liquid phase equilibrium variables and constraints
        ''' This section creates the necessary variables and constraints for
            calculating the liquid phase fugacity coefficient.
        '''
        blk = self
        
        # Compressibility factor
        blk.Z_liq_e = Var(domain=Reals, initialize=2.5e-2,
                    doc='Equilibrium liquid phase compressibility factor')

        # EoS coefficients at equilibrium pressure
        blk.A_liq_e = Var(domain=Reals,
                         doc='Liquid phase A coefficient at Pe')
        blk.B_liq_e = Var(domain=Reals,
                         doc='Liquid phase B coefficient at Pe')

        # Departure function terms
        blk.delta_liq = Var(blk.comp, domain=Reals,
                        doc='Coefficient for fugactiy departure functions')

        # ---------------------------------------------------------------------
        # Constraints
        # Equilibrium state coefficients
        blk.eq_A_liq_e = Constraint(expr = blk.A_liq_e*(blk.R*blk.T)**2
                                            == blk.am_liq*blk.Pe)
        blk.eq_B_liq_e = Constraint(expr = blk.B_liq_e*blk.R*blk.T
                                            == blk.bm_liq*blk.Pe)

        def rule_delta_liq(b, i):
            return b.delta_liq[i]*b.am_liq == 2*sqrt(blk.a[i]) \
                                    * sum(b.x[j]*sqrt(blk.a[j]) \
                                        * (1-b.kappa[i][j]) for j in b.comp)
        blk.eq_delta_liq = Constraint(blk.comp, rule=rule_delta_liq)

        # Equilibrium state compressibility factor calls
        def rule_Z_liq_e(b):
            return blk.proc_Z_liq(EoS_key[blk.EoS_ID],
                                          blk.A_liq_e,
                                          blk.B_liq_e)
        blk.expr_Z_liq_e = Expression(rule=rule_Z_liq_e)
        blk.eq_Z_liq_e = Constraint(expr = blk.Z_liq_e == blk.expr_Z_liq_e)


        # Fugacity coefficients
        def rule_phi_liq(b, j):
            return log(b.phi_liq[j])*(b.B_liq_e*b.EoS_p) \
                        == b.b[j]/b.bm_liq*(b.Z_liq_e-1)*(b.B_liq_e*b.EoS_p) \
                            - log(b.Z_liq_e-b.B_liq_e)*(b.B_liq_e*b.EoS_p) \
                            + b.A_liq_e * (b.b[j]/b.bm_liq - b.delta_liq[j]) \
                                * log((2*b.Z_liq_e + b.B_liq_e \
                                    *(b.EoS_u + b.EoS_p)) \
                                / (2*b.Z_liq_e + b.B_liq_e \
                                    * (b.EoS_u - b.EoS_p)))
        blk.eq_phi_liq = Constraint(blk.comp, rule=rule_phi_liq)

    def _make_equil_vap(self):
        # Create vapour phase equilibrium variables and constraints
        ''' This section creates the necessary variables and constraints for
            calculating the vapour phase fugacity coefficients.
        '''
        blk = self

        # Variables
        # Compressibility factors
        blk.Z_vap_e = Var(domain=Reals, initialize=0.85,
                    doc='Equilibrium vapour phase compressibility factor')

        # EoS coefficients at equilibrium pressure
        blk.A_vap_e = Var(domain=Reals,
                         doc='Vapour phase A coefficient at Pe')
        blk.B_vap_e = Var(domain=Reals,
                         doc='Vapour phase B coefficient at Pe')

        # Departure function terms
        blk.delta_vap = Var(blk.comp, domain=Reals,
                        doc='Coefficient for fugactiy departure functions')

        # ---------------------------------------------------------------------
        # Constraints
        # Equilibrium state coefficients
        blk.eq_A_vap_e = Constraint(expr = blk.A_vap_e*(blk.R*blk.T)**2
                                            == blk.am_vap*blk.Pe)
        blk.eq_B_vap_e = Constraint(expr = blk.B_vap_e*blk.R*blk.T
                                            == blk.bm_vap*blk.Pe)

        def rule_delta_vap(b, i):
            return b.delta_vap[i]*b.am_vap == 2*sqrt(blk.a[i]) \
                                    * sum(b.y[j]*sqrt(blk.a[j]) \
                                        * (1-b.kappa[i][j]) for j in b.comp)
        blk.eq_delta_vap = Constraint(blk.comp, rule=rule_delta_vap)

        # Equilibrium state compressibility factor calls
        def rule_Z_vap_e(b):
            return blk.proc_Z_vap(EoS_key[blk.EoS_ID],
                                          blk.A_vap_e,
                                          blk.B_vap_e)
        blk.expr_Z_vap_e = Expression(rule=rule_Z_vap_e)
        blk.eq_Z_vap_e = Constraint(expr = blk.Z_vap_e == blk.expr_Z_vap_e)

        # Fugacity coefficients
        def rule_phi_vap(b, j):
            return log(b.phi_vap[j])*(b.B_vap_e*b.EoS_p) \
                        == b.b[j]/b.bm_vap*(b.Z_vap_e-1)*(b.B_vap_e*b.EoS_p) \
                            - log(b.Z_vap_e-b.B_vap_e)*(b.B_vap_e*b.EoS_p) \
                            + b.A_vap_e * (b.b[j]/b.bm_vap - b.delta_vap[j]) \
                                * log((2*b.Z_vap_e + b.B_vap_e \
                                    *(b.EoS_u + b.EoS_p)) \
                                / (2*b.Z_vap_e + b.B_vap_e \
                                    * (b.EoS_u - b.EoS_p)))
        blk.eq_phi_vap = Constraint(blk.comp, rule=rule_phi_vap)

    def _make_ideal_props(self):
        # Create constraints and variables for ideal gas properties
        ''' This section creates the necessary variables and constraints for
            calculating the ideal gas mixture properties.
        '''
        blk = self
        
        # Variables
        blk.cp_ig_pc = Var(blk.comp, domain=Reals,
                    doc='Ideal gas pure component heat capacities [J/mol.K]')
        blk.h_ig_pc = Var(blk.comp, domain=Reals,
                    doc='Ideal gas pure component enthalpies [J/mol]')
        blk.s_ig_pc = Var(blk.comp, domain=Reals,
                    doc='Ideal gas pure component entropies [J/mol.K]')

        if blk.flag_phases != 'V':
            blk.cp_ig_liq = Var(domain=Reals,
                    doc='Liquid phase ideal gas heat capacity[J/mol.K]')
            blk.h_ig_liq = Var(domain=Reals,
                    doc='Liquid phase ideal gas enthalpy [J/mol]')
            blk.s_ig_liq = Var(domain=Reals,
                    doc='Liquid phase ideal gas entropy [J/mol.K]')
        if blk.flag_phases != 'L':
            blk.cp_ig_vap = Var(domain=Reals,
                    doc='Vapour phase ideal gas heat capacity [J/mol.K]')
            blk.h_ig_vap = Var(domain=Reals,
                    doc='Vapour phase ideal gas enthalpy [J/mol]')
            blk.s_ig_vap = Var(domain=Reals,
                    doc='Liquid phase ideal gas entropy [J/mol]')

        blk.P_vap = Var(blk.comp, domain=Reals, initialize=101325,
                    doc='Component ideal vapour pressure [Pa]')

        # ---------------------------------------------------------------------
        # Constraints
        # Ideal gas properties
        if blk.flag_phases != 'V':
            blk.eq_cp_ig_liq = Constraint(expr = blk.cp_ig_liq 
                            == sum(blk.x[j]*blk.cp_ig_pc[j] for j in blk.comp))
            blk.eq_h_ig_liq = Constraint(expr = blk.h_ig_liq 
                            == sum(blk.x[j]*blk.h_ig_pc[j] for j in blk.comp))
            blk.eq_s_ig_liq = Constraint(expr = blk.s_ig_liq 
                            == sum(blk.x[j]*blk.s_ig_pc[j] for j in blk.comp))

        if blk.flag_phases != 'L':
            blk.eq_cp_ig_vap = Constraint(expr = blk.cp_ig_vap
                            == sum(blk.y[j]*blk.cp_ig_pc[j] for j in blk.comp))
            blk.eq_h_ig_vap = Constraint(expr = blk.h_ig_vap
                            == sum(blk.y[j]*blk.h_ig_pc[j] for j in blk.comp))
            blk.eq_s_ig_vap = Constraint(expr = blk.s_ig_vap
                            == sum(blk.y[j]*blk.s_ig_pc[j] for j in blk.comp))

        def rule_P_vap(b, j):
            return log10(b.P_vap[j]*1e-5)*(b.T + b.Antoine[j][3]) \
                                == b.Antoine[j][1]*(b.T + b.Antoine[j][3]) \
                                    - b.Antoine[j][2]
        blk.eq_P_vap = Constraint(blk.comp, rule=rule_P_vap)

    def _make_liq_props(self):
        # Create thermodynamic variables and constraints
        ''' This section creates the necessary variables and constraints for
            calculating the thermodynamic properties of the liquid phase.
        '''
        blk = self

        # Variables
        # Additional pressure variables for smooth behaviour
        blk.Pl = Var(domain=Reals, initialize=101325,
                    doc='Liquid phase pressure [Pa]')

        # Phase compressibility factors
        blk.Z_liq = Var(domain=Reals, initialize=2.5e-2,
                    doc='Liquid phase compressibility factor')

        # Liquid phase properties
        blk.g_liq = Var(domain=Reals,
                        doc='Liquid phase molar Gibbs energy [J/mol]')
        blk.h_liq = Var(domain=Reals,
                        doc='Liquid phase molar enthalpy [J/mol]')
        blk.s_liq = Var(domain=Reals,
                        doc='Liquid phase molar entropy [J/mol.K]')
        blk.V_liq = Var(domain=Reals, initialize=1e-5,
                        doc='Liquid phase molar volume [m^3/mol]')

        # EoS coefficients at phase pressures
        blk.A_liq = Var(domain=Reals,
                         doc='Liquid phase A coefficient at Pl')
        blk.B_liq = Var(domain=Reals,
                         doc='Liquid phase B coefficient at Pl')

        # Departure function terms
        blk.dadT_liq = Var(domain=Reals, doc='Departure function coefficient')
        
        # ---------------------------------------------------------------------
        # Constraints
        # Equate inlet and phase mole fractions if one phase
        if blk.flag_phases != 'VL':
            def rule_zx(b, j):
                return b.x[j] == b.z[j]
            blk.eq_zx = Constraint(blk.comp, rule=rule_zx)

        # Phase state coefficients
        blk.eq_A_liq = Constraint(expr = blk.A_liq*(blk.R*blk.T)**2
                                            == blk.am_liq*blk.Pl)
        blk.eq_B_liq = Constraint(expr = blk.B_liq*blk.R*blk.T
                                            == blk.bm_liq*blk.Pl)

        blk.eq_dadT_liq = Constraint(expr = blk.dadT_liq*sqrt(blk.T) \
                        == -(blk.R/2)*sqrt(blk.omegaA) \
                        * sum(sum(blk.x[i]*blk.x[j]*(1-blk.kappa[i][j]) \
                        * (blk.fw[j]*sqrt(blk.a[i]*blk.Tc[j]/blk.Pc[j]) \
                        + blk.fw[i]*sqrt(blk.a[j]*blk.Tc[i]/blk.Pc[i])) \
                            for j in blk.comp) for i in blk.comp))

        # Phase pressure calculations
        if blk.flag_phases == 'VL':
            blk.eq_Pl = Constraint(expr = blk.Pl == blk.P \
                                + 0.5*(blk.Pe - blk.P \
                                    + sqrt((blk.Pe-blk.P)**2 + blk.eps**2)))
        else:
            blk.eq_Pl = Constraint(expr = blk.Pl == blk.P)

        # Phase compressibility factor calls
        def rule_Z_liq(b):
            return blk.proc_Z_liq(EoS_key[blk.EoS_ID],
                                          blk.A_liq,
                                          blk.B_liq)
        blk.expr_Z_liq = Expression(rule=rule_Z_liq)
        blk.eq_Z_liq = Constraint(expr = blk.Z_liq == blk.expr_Z_liq)

        # Liquid properties
        blk.eq_g_liq = Constraint(expr = blk.g_liq
                                        == blk.h_liq - blk.T*blk.s_liq)

        blk.eq_h_liq = Constraint(expr = (blk.h_liq - blk.h_ig_liq) \
                                * blk.bm_liq*blk.EoS_p \
                            == (blk.T*blk.dadT_liq - blk.am_liq) \
                            * log((2*blk.Z_liq + blk.B_liq \
                                    *(blk.EoS_u + blk.EoS_p)) \
                                / (2*blk.Z_liq + blk.B_liq \
                                    * (blk.EoS_u - blk.EoS_p))) \
                            + blk.R*blk.T*(blk.Z_liq-1)*(blk.bm_liq*blk.EoS_p))
        
        blk.eq_s_liq = Constraint(expr = (blk.s_liq - blk.s_ig_liq) \
                                * blk.bm_liq*blk.EoS_p \
                            == blk.R*log((blk.Z_liq-blk.B_liq)/blk.Z_liq) \
                                * blk.bm_liq*blk.EoS_p \
                            + blk.R*log(blk.Z_liq*blk.Pref/blk.Pl) \
                                * blk.bm_liq*blk.EoS_p \
                            + blk.dadT_liq * log((2*blk.Z_liq + blk.B_liq \
                                    *(blk.EoS_u + blk.EoS_p)) \
                                / (2*blk.Z_liq + blk.B_liq \
                                    * (blk.EoS_u - blk.EoS_p))))

        blk.eq_V_liq = Constraint(expr = blk.Z_liq*blk.T*blk.R
                                            == blk.Pl*blk.V_liq)

    def _make_vap_props(self):
        # Create thermodynamic variables and constraints
        ''' This section creates the necessary variables and constraints for
            calculating the thermodynamic properties of the vapour phase.
        '''
        blk = self

        # Variables
        # Additional pressure variables for smooth behaviour
        blk.Pv = Var(domain=Reals, initialize=101325,
                    doc='Vapour phase pressure [Pa]')

        # Phase compressibility factors
        blk.Z_vap = Var(domain=Reals, initialize=0.85,
                    doc='Vapour phase compressibility factor')

        # Vapour phase properties
        blk.g_vap = Var(domain=Reals,
                    doc='Vapour phase molar Gibbs Energy [J/mol]')
        blk.h_vap = Var(domain=Reals,
                    doc='Vapour phase molar enthalpy [J/mol]')
        blk.s_vap = Var(domain=Reals,
                    doc='Vapour phase molar entropy [J/mol.K]')
        blk.V_vap = Var(domain=Reals, initialize=1e-3,
                    doc='Vapour phase molar volume [m^3/mol]')

        # EoS coefficients at phase pressures
        blk.A_vap = Var(domain=Reals,
                         doc='Vapour phase A coefficient at Pv')
        blk.B_vap = Var(domain=Reals,
                         doc='Vapour phase B coefficient at Pv')

        # Departure function terms
        blk.dadT_vap = Var(domain=Reals, doc='Departure function coefficient')

        # ---------------------------------------------------------------------
        # Constraints
        # Equate inlet and phase mole fractions if one phase
        if blk.flag_phases != 'VL':
            def rule_zy(b, j):
                return b.y[j] == b.z[j]
            blk.eq_zy = Constraint(blk.comp, rule=rule_zy)

        # Phase state coefficients
        blk.eq_A_vap = Constraint(expr = blk.A_vap*(blk.R*blk.T)**2
                                            == blk.am_vap*blk.Pv)
        blk.eq_B_vap = Constraint(expr = blk.B_vap*blk.R*blk.T
                                            == blk.bm_vap*blk.Pv)

        blk.eq_dadT_vap = Constraint(expr = blk.dadT_vap*sqrt(blk.T)
                        == -(blk.R/2)*sqrt(blk.omegaA) \
                        * sum(sum(blk.y[i]*blk.y[j]*(1-blk.kappa[i][j]) \
                        * (blk.fw[j]*sqrt(blk.a[i]*blk.Tc[j]/blk.Pc[j]) \
                        + blk.fw[i]*sqrt(blk.a[j]*blk.Tc[i]/blk.Pc[i])) \
                            for j in blk.comp) for i in blk.comp))

        # Phase pressure calculations
        if blk.flag_phases == 'VL':
            blk.eq_Pv = Constraint(expr = blk.Pv == blk.P \
                                - 0.5*(blk.P - blk.Pe \
                                    + sqrt((blk.P-blk.Pe)**2 + blk.eps**2)))
        else:
            blk.eq_Pv = Constraint(expr = blk.Pv == blk.P)

        # Phase compressibility factor calls
        def rule_Z_vap(b):
            return blk.proc_Z_vap(EoS_key[blk.EoS_ID],
                                          blk.A_vap,
                                          blk.B_vap)
        blk.expr_Z_vap = Expression(rule=rule_Z_vap)
        blk.eq_Z_vap = Constraint(expr = blk.Z_vap == blk.expr_Z_vap)

        # Vapour properties
        blk.eq_g_vap = Constraint(expr = blk.g_vap
                                        == blk.h_vap - blk.T*blk.s_vap)

        blk.eq_h_vap = Constraint(expr = (blk.h_vap - blk.h_ig_vap) \
                                * blk.bm_vap*blk.EoS_p \
                            == (blk.T*blk.dadT_vap - blk.am_vap) \
                            * log((2*blk.Z_vap + blk.B_vap \
                                    *(blk.EoS_u + blk.EoS_p)) \
                                / (2*blk.Z_vap + blk.B_vap \
                                    * (blk.EoS_u - blk.EoS_p))) \
                            + blk.R*blk.T*(blk.Z_vap-1)*blk.bm_vap*blk.EoS_p)

        blk.eq_s_vap = Constraint(expr = (blk.s_vap - blk.s_ig_vap) \
                                * blk.bm_vap*blk.EoS_p\
                            == blk.R*log((blk.Z_vap-blk.B_vap)/blk.Z_vap) \
                                * blk.bm_vap*blk.EoS_p \
                            + blk.R*log(blk.Z_vap*blk.Pref/blk.Pv) \
                                * blk.bm_vap*blk.EoS_p \
                            + blk.dadT_vap * log((2*blk.Z_vap + blk.B_vap \
                                    *(blk.EoS_u + blk.EoS_p)) \
                                / (2*blk.Z_vap + blk.B_vap \
                                    * (blk.EoS_u - blk.EoS_p))))

        blk.eq_V_vap = Constraint(expr = blk.Z_vap*blk.T*blk.R
                                            == blk.Pv*blk.V_vap)

    def _make_thermo_mix(self):
        # Create variables and constraints for mixture properties
        ''' This section creates the necessary variables and constraints for
            calculating the thermodynamic properties in the mixture.
        '''
        blk = self
        
        # Variables
        blk.g_mix = Var(domain=Reals, doc='Mixture molar Gibbs Energy[J/mol]')
        blk.h_mix = Var(domain=Reals, doc='Mixture molar enthalpy [J/mol]')
        blk.s_mix = Var(domain=Reals, doc='Mixture molar entrolpy [J/mol.K]')
        blk.V_mix = Var(domain=Reals, doc='Mixture molar volume [m^3/mol]')

        # ---------------------------------------------------------------------
        # Constraints
        blk.eq_g_mix = Constraint(expr = blk.g_mix
                                        == blk.V*blk.g_vap + blk.L*blk.g_liq)
        blk.eq_h_mix = Constraint(expr = blk.h_mix
                                        == blk.V*blk.h_vap + blk.L*blk.h_liq)
        blk.eq_s_mix = Constraint(expr = blk.s_mix
                                        == blk.V*blk.s_vap + blk.L*blk.s_liq)
        blk.eq_V_mix = Constraint(expr = blk.V_mix
                                        == blk.V*blk.V_vap + blk.L*blk.V_liq)

    def _make_trans_liq(self):
        # Create constraints for liquid phase transport properties
        ''' This section creates the necessary variables and constraints
            for calculating the liquid phase transport properties.
        '''
        blk = self
        
        # Variables
        blk.mw_liq = Var(domain=Reals, initialize=0.03,
                         doc='Liquid phase mixture molecular weight [mol/kg]')

        # ---------------------------------------------------------------------
        # Constraints
        # Molecular weight
        blk.eq_mw_liq = Constraint(expr = blk.mw_liq
                            == sum(blk.x[j]*blk.mw_pc[j] for j in blk.comp))

    def _make_trans_vap(self):
        # Create constraints for vapour phase transport properties
        ''' This section creates the necessary variables and constraints
            for calculating the vapour phase transport properties.
        '''
        blk = self
        
        # Variables
        blk.mw_vap = Var(domain=Reals, initialize=0.03,
                         doc='Vapour phase mixture molecular weight [mol/kg]')

        # ---------------------------------------------------------------------
        # Constraints
        # Molecular weight
        blk.eq_mw_vap = Constraint(expr = blk.mw_vap
                            == sum(blk.y[j]*blk.mw_pc[j] for j in blk.comp))

    def _make_trans_mix(self):
        # Create constraints for mixture transport properties
        ''' This section creates the necessary variables and constraints
            for calculating the mixture transport properties.
        '''
        blk = self
        
        # Variables
        blk.mw_mix = Var(domain=Reals, initialize=0.03,
                         doc='Total mixture molecular weight [mol/kg]')

        # ---------------------------------------------------------------------
        # Constraints
        # Molecular weight
        blk.eq_mw_mix = Constraint(expr = blk.mw_mix
                            == sum(blk.z[j]*blk.mw_pc[j] for j in blk.comp))

    def _initialize(blk, T=None, P=None, x=None, y=None, z=None,
                        outlvl=0, optarg=None):
        ''' Initialisation routine for property package (default solver ipopt)
        
            Arguments:
                T - T value to use for initialisation (not used if T is fixed)
                P - P value to use for initialisation (not used if T is fixed)
                z - array of z values to use for initialisation (sum(z)=1)
                        (not used if z is fixed)
                x,y - arrays of initial guesses for phase compositions
                outlvl - sets output level of initialisation routine
                            0 = no output
                            1 = return solver state for each step in routine
                            2 = include solver output infomation (tee=True)
                optarg - solver options dictionary object
        '''
        # Fix state variables if not already fixed
        if blk.P.fixed == True:
            Pflag = 1
        else:
            Pflag = 0
            if P == None:
                blk.P.fix(101325)
            else:
                blk.P.fix(P)

        if blk.T.fixed == True:
            Tflag = 1
        else:
            Tflag = 0        
            if T == None:
                blk.T.fix(298.15)
            else:
                blk.T.fix(T)

        zflag={}
        for j in blk.comp:
            if blk.z[j].fixed == True:
                zflag[j] = 1
            else:
                zflag[j] = 0
                if z == None:
                    blk.z[j].fix(1/len(blk.comp))
                else:
                    blk.z[j].fix(z[j])

        # Set inital guesses for phase compositions
        if blk.flag_phases != 'V':
            if x == None:
                for j in blk.comp:
                    blk.x[j].fix(value(blk.z[j]))
            else:
                for j in blk.comp:
                    blk.x[j].fix(x[j])

        if blk.flag_phases != 'L':
            if y == None:
                for j in blk.comp:
                    blk.y[j].fix(value(blk.z[j]))
            else:
                for j in blk.comp:
                    blk.y[j].fix(y[j])

        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False
        
        if optarg == None:
            sopt = {'tol':1e-10}
        else:
            sopt=optarg

        # ---------------------------------------------------------------------
        # Deactivate constraints
        if blk.flag_phases == 'VL':
            blk.eq_Keq.deactivate()
            blk.eq_mbal_t.deactivate()
            blk.eq_sum_xy.deactivate()
            blk.eq_mbal_c.deactivate()
            blk.eq_Pe.deactivate()
            blk.eq_beta.deactivate()
            blk.eq_GL.deactivate()
            blk.eq_GV.deactivate()

            blk.eq_Z_liq_e.deactivate()
            blk.eq_delta_liq.deactivate()
            blk.eq_phi_liq.deactivate()

            blk.eq_Z_vap_e.deactivate()
            blk.eq_delta_vap.deactivate()
            blk.eq_phi_vap.deactivate()

        if blk.flag_phases != 'V':
            blk.eq_Pl.deactivate()
            blk.eq_Z_liq.deactivate()
            blk.eq_g_liq.deactivate()
            blk.eq_h_liq.deactivate()
            blk.eq_s_liq.deactivate()
            blk.eq_V_liq.deactivate()
            blk.eq_A_liq.deactivate()
            blk.eq_B_liq.deactivate()
            blk.eq_dadT_liq.deactivate()
            blk.eq_cp_ig_liq.deactivate()
            blk.eq_h_ig_liq.deactivate()
            blk.eq_s_ig_liq.deactivate()
            if blk.flag_phases != 'VL':
                blk.eq_zx.deactivate()
        
        if blk.flag_phases != 'L':
            blk.eq_Pv.deactivate()
            blk.eq_Z_vap.deactivate()
            blk.eq_g_vap.deactivate()
            blk.eq_h_vap.deactivate()
            blk.eq_s_vap.deactivate()
            blk.eq_V_vap.deactivate()
            blk.eq_A_vap.deactivate()
            blk.eq_B_vap.deactivate()
            blk.eq_dadT_vap.deactivate()
            blk.eq_cp_ig_vap.deactivate()
            blk.eq_h_ig_vap.deactivate()
            blk.eq_s_ig_vap.deactivate()
            if blk.flag_phases != 'VL':
                blk.eq_zy.deactivate()

        if blk.flag_mprop == 1:
            blk.eq_g_mix.deactivate()
            blk.eq_h_mix.deactivate()
            blk.eq_s_mix.deactivate()
            blk.eq_V_mix.deactivate()

        if blk.flag_phases != 'V':
            blk.eq_mw_liq.deactivate()
        
        if blk.flag_phases != 'L':
            blk.eq_mw_vap.deactivate()

        if blk.flag_mprop == 1:
            blk.eq_mw_mix.deactivate()

        # ---------------------------------------------------------------------
        # Step 1 - Initialise EoS variables
        if blk.flag_phases == 'VL':
            blk.Pe.fix(value(blk.P))

        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            if results.solver.status == SolverStatus.ok:
                print(blk,"Initialisation Step  1 Complete")
            else:
                print(blk,"Initialisation Step  1 Failed")

        # ---------------------------------------------------------------------
        # Step 2 - Initialise phase equilibrium 1
        if blk.flag_phases == 'VL':
            blk.Pe = sum(value(blk.z[j])*value(blk.P_vap[j]) for j in blk.comp)
            for i in blk.comp:
                blk.delta_liq[i] = (2*sqrt(value(blk.a[i]))/value(blk.am_liq))\
                            * sum(value(blk.x[j])*sqrt(value(blk.a[j])) \
                                * (1-value(blk.kappa[i][j])) for j in blk.comp)
                blk.delta_vap[i] = (2*sqrt(value(blk.a[i]))/value(blk.am_vap))\
                            * sum(value(blk.y[j])*sqrt(value(blk.a[j])) \
                                * (1-value(blk.kappa[i][j])) for j in blk.comp)

            blk.eq_delta_liq.activate()
            blk.eq_delta_vap.activate()

            results = blk.solve(options=sopt,tee=stee)
            if outlvl > 0:
                if results.solver.status == SolverStatus.ok:
                    print(blk,"Initialisation Step  2 Complete")
                else:
                    print(blk,"Initialisation Step  2 Failed")
        else:
            if outlvl > 0:
                print(blk,"Initialisation Step  2 Skipped")

        # ---------------------------------------------------------------------
        # Step 3 - Initialise phase equilibrium 2
        if blk.flag_phases == 'VL':
            blk.Z_liq_e = value(blk.expr_Z_liq_e)
            blk.Z_vap_e = value(blk.expr_Z_vap_e)

            blk.eq_phi_liq.activate()
            blk.eq_phi_vap.activate()

            results = blk.solve(options=sopt,tee=stee)
            if outlvl > 0:
                if results.solver.status == SolverStatus.ok:
                    print(blk,"Initialisation Step  3 Complete")
                else:
                    print(blk,"Initialisation Step  3 Failed")
        else:
            if outlvl > 0:
                print(blk,"Initialisation Step  3 Skipped")

        # ---------------------------------------------------------------------
        # Step 4 - Initialise phase equilibrium 2
        if blk.flag_phases == 'VL':
            blk.beta = value(blk.P)/value(blk.Pe)
            for j in blk.comp:
                blk.y[j] = (value(blk.z[j])*value(blk.beta) \
                            * value(blk.phi_liq[j]) / value(blk.phi_vap[j])) \
                / sum(value(blk.z[k])*value(blk.beta)*value(blk.phi_liq[k]) \
                            / value(blk.phi_vap[k]) for k in blk.comp)
                blk.x[j] = (value(blk.z[j])/(value(blk.beta) \
                            * value(blk.phi_liq[j]) / value(blk.phi_vap[j]))) \
                / sum(value(blk.z[k])/(value(blk.beta)*value(blk.phi_liq[k]) \
                            / value(blk.phi_vap[k])) for k in blk.comp)

            blk.Pe.fixed = False
            for i in blk.comp:
                blk.x[i].fixed = False
                blk.y[i].fixed = False

            blk.eq_Pe.activate()
            blk.eq_beta.activate()
            blk.eq_GL.activate()
            blk.eq_GV.activate()

            blk.eq_Z_liq_e.activate()
            blk.eq_Z_vap_e.activate()
        
            blk.eq_Keq.activate()
            blk.eq_mbal_t.activate()
            blk.eq_sum_xy.activate()
            blk.eq_mbal_c.activate()

            results = blk.solve(options=sopt,tee=stee)
            if outlvl > 0:
                if results.solver.status == SolverStatus.ok:
                    print(blk,"Initialisation Step  4 Complete")
                else:
                    print(blk,"Initialisation Step  4 Failed")
        else:
            if outlvl > 0:
                print(blk,"Initialisation Step  4 Skipped")

        # ---------------------------------------------------------------------
        # Step 5 - Initialise phase thermo properties 1
        if blk.flag_phases != 'V':
            if blk.flag_phases == 'VL':
                blk.Pl = max(value(blk.P),value(blk.Pe))
            else:
                blk.Pl = value(blk.P)

                blk.eq_zx.activate()
                for j in blk.comp:
                    blk.x[j].fixed = False

            blk.eq_A_liq.activate()
            blk.eq_B_liq.activate()
            blk.eq_dadT_liq.activate()
            blk.eq_Pl.activate()

            blk.eq_cp_ig_liq.activate()
            blk.eq_h_ig_liq.activate()
            blk.eq_s_ig_liq.activate()

        if blk.flag_phases != 'L':
            if blk.flag_phases == 'VL':
                blk.Pv = min(value(blk.P),value(blk.Pe))
            else:
                blk.Pv = value(blk.P)

                blk.eq_zy.activate()
                for j in blk.comp:
                    blk.y[j].fixed = False

            blk.eq_A_vap.activate()
            blk.eq_B_vap.activate()
            blk.eq_dadT_vap.activate()
            blk.eq_Pv.activate()

            blk.eq_cp_ig_vap.activate()
            blk.eq_h_ig_vap.activate()
            blk.eq_s_ig_vap.activate()

        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            if results.solver.status == SolverStatus.ok:
                print(blk,"Initialisation Step  5 Complete")
            else:
                print(blk,"Initialisation Step  5 Failed")

        # ---------------------------------------------------------------------
        # Step 6 - Initialise phase thermo properties 2
        if blk.flag_phases != 'V':
            blk.Z_liq = value(blk.expr_Z_liq)

            blk.eq_Z_liq.activate()
            blk.eq_g_liq.activate()
            blk.eq_h_liq.activate()
            blk.eq_s_liq.activate()
            blk.eq_V_liq.activate()

        if blk.flag_phases != 'L':
            blk.Z_vap = value(blk.expr_Z_vap)

            blk.eq_Z_vap.activate()
            blk.eq_g_vap.activate()
            blk.eq_h_vap.activate()
            blk.eq_s_vap.activate()
            blk.eq_V_vap.activate()
        
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            if results.solver.status == SolverStatus.ok:
                print(blk,"Initialisation Step  6 Complete")
            else:
                print(blk,"Initialisation Step  6 Failed")

        # ---------------------------------------------------------------------
        # Step 7 - Initialise mixture thermo properties
        if blk.flag_mprop == 1:
            blk.eq_g_mix.activate()
            blk.eq_h_mix.activate()
            blk.eq_s_mix.activate()
            blk.eq_V_mix.activate()

            results = blk.solve(options=sopt,tee=stee)
            if outlvl > 0:
                if results.solver.status == SolverStatus.ok:
                    print(blk,"Initialisation Step  7 Complete")
                else:
                    print(blk,"Initialisation Step  7 Failed")
        else:
            if outlvl > 0:
                print(blk,"Initialisation Step  7 Skipped")

        # ---------------------------------------------------------------------
        # Step 8 - Initialise phase transport properties 1
        if blk.flag_phases != 'V':
            blk.eq_mw_liq.activate()

        if blk.flag_phases != 'L':
            blk.eq_mw_vap.activate()
        
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            if results.solver.status == SolverStatus.ok:
                print(blk,"Initialisation Step  8 Complete")
            else:
                print(blk,"Initialisation Step  8 Failed")

        # ---------------------------------------------------------------------
        # Step 9 - Initialise mixture transport properties
        if blk.flag_mprop == 1:
            blk.eq_mw_mix.activate()

            results = blk.solve(options=sopt,tee=stee)
            if outlvl > 0:
                if results.solver.status == SolverStatus.ok:
                    print(blk,"Initialisation Step  9 Complete")
                else:
                    print(blk,"Initialisation Step  9 Failed")
        else:
            if outlvl > 0:
                print(blk,"Initialisation Step  9 Skipped")

        # ---------------------------------------------------------------------
        # Unfix state variables
        if Pflag == 0:
            blk.P.fixed = False

        if Tflag == 0:
            blk.T.fixed = False

        for j in blk.comp:
            if zflag[j] == 0:
                blk.z[j].fixed = False

        if outlvl > 0:
            print(blk,"Initialisation Complete")

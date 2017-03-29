"""
Bubbling Fluidised Bed Reactor Model

Add documentation here

To Do:
Collocation discretisation
fc_max limit - code fragment present in comments
Switch to two step bulk gas transfer
Bubble diamter constraint
Slugging
Convert solids composition from loading to mass fraction
Solids momentum balances
Reversible gas flow in emulsion region

"""
from __future__ import division
from __future__ import print_function
__author__ = "Andrew Lee"
__version__ = "7.0.0 r1"

from pyomo.environ import *
from pyomo.dae import *

import time

# Import IDAES cores
from idaes_models.core import UnitModel, ProcBlock
import idaes_models.core.util.misc as model_util

__all__ = ['BFB Reactor']

@ProcBlock("BFB")
class _BFB(UnitModel):
    def __init__(self, *args, **kwargs):
        """
        Args:
        fe_set = set of normalised finite element locations
        nfe = number of finite elements for bed discretization
                    (not used if fe_set specified)
        s_inlet = type of solids feed (Top or Bottom)
        s_outlet = type of solids outlet (Overflow or Underflow)
        hx_type = arrangement of heat exchanger tubes (not used yet)
        gas_prop_lib = library file for gas phase properties
        sol_prop_lib = library file for solid phase properties
        hx_prop_lib = library file for heat ecxhanger fluid properties
        vb_method = correlation to use to calculate bubble velocity
                    - Davidson
                    - Werther A
                    - Werther B
        """
        self.fe_set_temp = kwargs.pop("fe_set",None)
        
        if self.fe_set_temp == None:
            self.nfe = kwargs.pop("nfe", 60)
            self.fe_set = {}
            for i in range(0,self.nfe+1):
                self.fe_set[float(i)] = float(i/self.nfe)
        else:
            self.fe_set = self.fe_set_temp
            self.nfe = len(self.fe_set) - 1
            kwargs.pop("nfe", 60)

        self.s_inlet = kwargs.pop("s_inlet", "Bottom")
        if not((self.s_inlet == 'Top') or (self.s_inlet == 'Bottom')):
            raise Exception('Solid feed type not recognised. '\
                            'Must be either Top or Bottom')

        self.s_outlet = kwargs.pop("s_outlet", "Overflow") 
        if not ((self.s_outlet=='Overflow') or (self.s_outlet=='Underflow')):
            raise Exception('Solid outlet type not recognised. '\
                            'Must be either Overflow or Underflow')

        self.hx_type = kwargs.pop("hx_type", "Tubes")
        if not((self.hx_type == 'Tubes')):
            raise Exception('HX type not recognised. '\
                            'Must be Tubes (for now)')

        self.gas_prop_lib = kwargs.pop("gas_prop_lib")
        self.sol_prop_lib = kwargs.pop("sol_prop_lib")
        self.hx_prop_lib = kwargs.pop("hx_prop_lib")
        
        self.vb_method = kwargs.pop("vb_method","Davidson")

        UnitModel.__init__(self, *args, **kwargs)

    def build(self, *args, **kwargs):

        self._import_props()

        self._make_params()
        self._make_domain()

        self._make_vars()
        self._make_hydro_vars()
        self._make_HX_vars()

        self._dae_transform()

        self._make_props()

        self._make_bed_model()
        self._make_bdry_conds()
        self._make_hydro_model()
        self._make_HX_model()

        return self

    def _import_props(self):
        # Import property packages
        blk = self
    
        blk.gas_prop_lib = __import__(self.gas_prop_lib)
        blk.sol_prop_lib = __import__(self.sol_prop_lib)
        blk.hx_prop_lib = __import__(self.hx_prop_lib)

        # Declare Component Lists
        blk.GasList = blk.gas_prop_lib.comp
        blk.SolidList = blk.sol_prop_lib.comp
        blk.HXList = blk.hx_prop_lib.comp 
        
    def _make_params(self):
        """
        Create the paramters and sets for the Pyomo block
        """
        blk = self

        # Declare Imutable Parameters
        blk.pi = Param(default=2*acos(0), doc='pi')
        blk.R  = Param(default=8.314472,
                      doc='Universal Gas Constant [J/mol.K]')
        blk.gc = Param(default=9.81,
                      doc='Gravitational Acceleration Constant [m^2/s]')
        blk.eps = Param(default=1e-4,
                      doc='Smoothing Factor for Smooth IF Statements')

        # Declare Mutable Parameters
        blk.ah = Param(within=PositiveReals, mutable=True, default=0.8,
                      doc='Emprical Factor for HX Tube Model')
        blk.Cr = Param(within=PositiveReals, mutable=True, default=1,
                      doc='Average Correction Factor for HX Tube Model')
        blk.fc_max = Param(within=PositiveReals, mutable=True, default=0.1,
                      doc='Maximum Bubble to Cloud Volume Ratio')
        blk.fw = Param(within=PositiveReals, mutable=True, default=0.2,
                      doc='Bubble to Wake Volume Ratio')
        blk.hw = Param(within=PositiveReals, mutable=True, default=1500,
                    doc='Heat Transfer Coefficient of HX Tube Walls [W/m^2.K]')
        blk.Kd = Param(within=NonNegativeReals, mutable=True, default=1,
                      doc='Bulk Gas Permeation Coefficient [m/s]')

    def _make_domain(self):
        """
        Create the axial domains for the Pyomo block
        """
        blk = self

        # Declare Distribution Domain
        blk.l = ContinuousSet(bounds=(0.0,float(blk.nfe)),
                    doc='Set of Discrete Elements')

    def _make_vars(self):
        """
        Create the variables for the Pyomo block
        """
        blk = self

        # Discterisation point location and sizeSize
        blk.ll = Var(blk.l, domain=Reals,
                     doc='Location of Discrete Elements [m]')   
        
        # Vessel Dimensions
        blk.Ax = Var(domain=Reals,
                    doc='Cross-Sectional Area of Bed [m^2]')
        blk.Areact = Var(domain=Reals,
                    doc='Cross-Sectional Area of Reactor Vessel [m^2]')
        blk.Dt = Var(domain=Reals,
                    doc='Diameter of Reactor Vessel [m]')
        blk.Dte = Var(domain=Reals,
                    doc='Hydraulic Diameter of Bed [m]')
        blk.Lb = Var(domain=Reals,
                    doc='Depth of Bed [m]')

        # Distributor Design
        blk.Ao = Var(domain=Reals,
                    doc='Distributor Plate Area per Orifice [m^2/orifice]')
        blk.nor = Var(domain=Reals,
                      doc='Distributor Plate Orifices per Area [orifices/m^2]')

        # Gas Inlet Conditions
        blk.Gas_In_F = Var(domain=Reals,
                    doc='Molar Flowrate of Gas at Gas Inlet [mol/s]')
        blk.Gas_In_P = Var(domain=Reals,
                    doc='Pressure of Gas at Gas Inlet [Pa]')
        blk.Gas_In_T = Var(domain=Reals,
                    doc='Temperature of Gas at Gas Inlet [K]')
        blk.Gas_In_y = Var(blk.GasList, domain=Reals,
                    doc='Mole Fractions of Gas Species at Gas Inlet [mol/mol]')

        # Solids Inlet Conditions
        blk.Solid_In_F = Var(domain=Reals,
                    doc='Mass Flowrate of Solids at Solid Inlet [kg/s]')
        blk.Solid_In_T = Var(domain=Reals,
                    doc='Temperature of Solids at Solids Inlet [k]')
        blk.Solid_In_x = Var(blk.SolidList, domain=Reals,
                    doc='Active Mass of Solids at Solids Inlet [kg/kg]')

        # Material Flows
        blk.Gb = Var(blk.l, domain=Reals,
                    doc='Bubble Region Molar Gas FLowrate [mol/s]')
        blk.Ge = Var(blk.l, domain=Reals,
                    doc='Emulsion Region Molar Gas FLowrate [mol/s]')
        blk.Jc = Var(blk.l, domain=Reals,
                    doc='Cloud-Wake Region Solids Flux [kg/m^2.s]')
        blk.Je = Var(blk.l, domain=Reals,
                    doc='Emulsion Region Solids Flux (downwards) [kg/m^2.s]')

        # Component Material Flows
        blk.Gbc = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Bubble Region Component Molar Gas FLowrate [mol/s]')
        blk.Gec = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Emulsion Region Component Molar Gas FLowrate [mol/s]')
        blk.Jcc = Var(blk.SolidList, blk.l, domain=Reals,
                    doc='Cloud-Wake Region Adsorbed Species Flux [mol/m^2.s]')
        blk.Jec = Var(blk.SolidList, blk.l, domain=Reals,
                    doc='Emulsion Region Adsorbed Species Flux [mol/m^2.s]')

        # Enthalpy Flows
        blk.Gbh = Var(blk.l,doc='Bubble Region Gas Enthalpy Flowrate [J/s]')
        blk.Geh = Var(blk.l,doc='Emulsion Region Gas Enthalpy Flowrate [J/s]')
        blk.Jch = Var(blk.l,
                    doc='Cloud-Wake Region Solids Enthalpy Flux [J/m^2.s]')
        blk.Jeh = Var(blk.l,
                    doc='Emulsion Region Solids Enthalpy Flux [J/m^2.s]')

        # Temperatures and Pressures
        blk.P = Var(blk.l, domain=Reals, doc='Pressure [Pa]')
        blk.Tgb = Var(blk.l, domain=Reals,
                      doc='Bubble Region Gas Temperature [K]')
        blk.Tgc = Var(blk.l, domain=Reals,
                      doc='Cloud-Wake Region Gas Temperature [K]')
        blk.Tge = Var(blk.l, domain=Reals,
                      doc='Emulsion Region Gas Temperature [K]')
        blk.Tsc = Var(blk.l, domain=Reals,
                      doc='Cloud-Wake Region Solids Temperature [K]')
        blk.Tse = Var(blk.l, domain=Reals,
                      doc='Emulsion Region Solids Temperature [K]')

        # Gas and Solid Compositions
        blk.yb = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Bubble Region Gas Mole Fractions [mol/mol]')
        blk.yc = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Cloud-Wake Region Gas Mole Fractions [mol/mol]')
        blk.ye = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Emulsion Region Gas Mole Fractions [mol/mol]')
        blk.xc = Var(blk.SolidList, blk.l, domain=Reals,
                    doc='Cloud-Wake Region Solids Loading [mol/kg]')
        blk.xe = Var(blk.SolidList, blk.l, domain=Reals,
                    doc='Emulsion Region Solids Loading [mol/kg')

        # Gas Phase Concentrations
        blk.cb = Var(blk.GasList, blk.l, domain=Reals,
                doc='Bubble Region Gas Component Concentrations [mol/m^3]')
        blk.cc = Var(blk.GasList, blk.l, domain=Reals,
                doc='Cloud-Wake Region Gas Component Concentrations [mol/m^3]')
        blk.ce = Var(blk.GasList, blk.l, domain=Reals,
                doc='Emulsion Region Gas Component Concentrations [mol/m^3]')
        blk.cct = Var(blk.l, domain=Reals,
                doc='Cloud-Wake Region Total Gas Concentration [mol/m^3]')
        blk.cet = Var(blk.l, domain=Reals,
                doc='Emulsion Region Total Gas Concentration [mol/m^3]')

        # Velocities
        blk.vb = Var(blk.l, domain=Reals,
                    doc='Bubble Velocity [m/s]')
        blk.ve = Var(blk.l, doc='Emulsion Region Gas Velocity [m/s]')
        blk.vg = Var(blk.l, domain=Reals,
                    doc='Superficial Gas Veleocity [m/s]')
        blk.us = Var(blk.l,
                    doc='Emulsion Region Solids Velocity (downwards) [m/s]')

        # Fluidisation Properties
        blk.Ar = Var(blk.l, domain=Reals,
                    doc='Archimedes Number [-]')

        # Bubble Dimensions and Hydrodynamics
        blk.db = Var(blk.l, domain=Reals,
                    doc='Average Bubble Diameter [m]')
        blk.delta = Var(blk.l, domain=Reals,
                    doc='Bubble Region Volume Fraction [m^3/m^3]')
        blk.delta_e = Var(blk.l, domain=Reals,
                    doc='Emulsion Region Volume Fraction [m^3/m^3]')
        blk.e = Var(blk.l, domain=Reals,
                    doc='Cross-Sectional Average Voidage Fraction [m^3/m^3]')
        blk.ed = Var(blk.l, domain=Reals,
                    doc='Emulsion Region Voidage Fraction [m^3/m^3]')
        blk.fc = Var(blk.l, domain=Reals,
                    doc='Cloud to Bubble Region Volume Ratio [m^3/m^3]')
        blk.fcw = Var(blk.l, domain=Reals,
                    doc='Cloud-Wake to Bubble Region Volume Ratio [m^3/m^3]')

        # Bulk Material Transfer
        blk.Kgbulk_c = Var(blk.GasList, blk.l,
                    doc='Gas Phase Component Bulk Transfer Rate [mol/s]')
        blk.Ksbulk = Var(blk.l, doc='Solid Phase Bulk Transfer Rate [kg/s]')
        blk.Ksbulk_c = Var(blk.SolidList, blk.l,
                    doc='Adsorbed Species Bulk Transfer Rate [mol/s]')
        blk.Hgbulk = Var(blk.l,
                    doc='Gas Phase Bulk Enthalpy Trasnfer Rate [J/s]')
        blk.Hsbulk = Var(blk.l,
                    doc='Solid Phase Bulk Enthalpy Trasnfer Rate [J/s]')

        # Heat and Mass Transfer Coefficients
        blk.Kbc = Var(blk.GasList, blk.l, domain=Reals,
                    doc='Bubble to Cloud-Wake Gas Mass Transfer Coefficient '\
                    '[1/s]')
        blk.Kce = Var(blk.GasList, blk.l, domain=Reals,
                     doc='Cloud-Wake to Emulsion Gas Mass Transfer '\
                     'Coefficient [1/s]')
        blk.Kcebs = Var(blk.l,
                    doc='Cloud-Wake to Emulsion Solid Mass Transfer '\
                    'Coefficient [1/s]')
        blk.Kcebs_c = Var(blk.SolidList, blk.l, domain=Reals,
                    doc='Cloud-Wake to Emulsion Adsorbed Species Mass '\
                    'Transfer Rate [mol/s]')
        blk.Hbc = Var(blk.l, domain=Reals,
                    doc='Bubble to Cloud-Wake Gas Energy Transfer '\
                    'Coefficient [J/m^3.K.s]')
        blk.Hce = Var(blk.l, domain=Reals,
                    doc='Cloud-Wake to Emulsion Gas Energy Transfer '\
                    'Coefficient [J/m^3.K.s]')
        blk.Hcebs = Var(blk.l,
                    doc='Cloud-Wake to Emulsion Solid Enthalpy Transfer '\
                    'Rate [J/s]')
        blk.hp = Var(blk.l, domain=Reals,
                    doc='Gas to Solid Energy Convective Heat Transfer '\
                    'Coefficient [J/m^2.K.s]')
        blk.hp_c = Var(blk.l, domain=Reals,
                    doc='Cloud-Wake Region Gas to Solid Enthalpy Transfer '\
                    'Rate [J/s]')
        blk.hp_e = Var(blk.l, domain=Reals,
                    doc='Emulsion Region Gas to Solid Enthalpy Transfer '\
                    'Rate [J/s]')
        blk.hr_c = Var(blk.l, domain=Reals,
                       doc='Cloud-Wake Region Gas-Solids Heat Transfer due '\
                       'to Reaction [J/s]')
        blk.hr_e = Var(blk.l, domain=Reals,
                       doc='Emulsion Region Gas-Solids Heat Transfer due '\
                       'to Reaction [J/s]')

        # Reaction Rates
        blk.rgc = Var(blk.GasList, blk.l,
                      doc='Cloud-Wake Region Gas Phase Reaction Rates [mol/s]')
        blk.rge = Var(blk.GasList, blk.l,
                      doc='Emulsion Region Gas Phase Reaction Rates [mol/s]')
        blk.rsc = Var(blk.SolidList, blk.l,
                      doc='Cloud-Wake Region Solid Phase Reaction Rates '\
                      '[mol/s]')
        blk.rse = Var(blk.SolidList, blk.l,
                      doc='Cloud-Wake Region Solid Phase Reaction Rates '\
                      '[mol/s]')

        # Heat Transfer from Bed
        blk.Q = Var(doc='Total Heat Transfered from Bed [J/s]')
        blk.Qhx = Var(blk.l,
                    doc='Heat Transfered to from Bed in Discrete Element '\
                    '[J/s]')

        # Outlet Variables
        # Gas Outlet Conditions
        blk.Gas_Out_F = Var(domain=Reals,
                    doc='Gas Flowrate at Gas Outlet [mol/s]')
        blk.Gas_Out_P = Var(domain=Reals,
                    doc='Pressure at Gas Outlet [Pa]')
        blk.Gas_Out_T = Var(domain=Reals,
                    doc='Gas Temeprature at Gas Outlet [K]')
        blk.Gas_Out_y = Var(blk.GasList, domain=Reals,
                    doc='Gas Mole Fractions at Gas Outlet [mol/mol]')

        # Solids Inlet Conditions
        blk.Solid_Out_F = Var(domain=Reals,
                    doc='Solid Flowrate at Solid Outlet [kg/s]')
        blk.Solid_Out_T = Var(domain=Reals,
                    doc='Solid Temperature at Solid Outlet [K]')
        blk.Solid_Out_x = Var(blk.SolidList, domain=Reals,
                    doc='Solid Active Mass at Solid Outlet [kg/kg]')

        # Derivative Variables
        blk.dGbcdx = DerivativeVar(blk.Gbc, wrt=blk.l)
        blk.dGecdx = DerivativeVar(blk.Gec, wrt=blk.l)
        blk.dGbhdx = DerivativeVar(blk.Gbh, wrt=blk.l)
        blk.dGehdx = DerivativeVar(blk.Geh, wrt=blk.l)
        blk.dJcdx = DerivativeVar(blk.Jc, wrt=blk.l)
        blk.dJccdx = DerivativeVar(blk.Jcc, wrt=blk.l)
        blk.dJecdx = DerivativeVar(blk.Jec, wrt=blk.l)
        blk.dJchdx = DerivativeVar(blk.Jch, wrt=blk.l)
        blk.dJehdx = DerivativeVar(blk.Jeh, wrt=blk.l)
        blk.dPdx = DerivativeVar(blk.P, wrt=blk.l)
        blk.dl = DerivativeVar(blk.ll, wrt=blk.l,
                    doc='Size of Discrete Elements [m]')

    def _make_bed_model(self):
        """
        Create main body of model.
        
        Groups: A,B (except boundaries), C, D, E, F, G, L
        """
        blk = self

        # Reactor Sizing and Design
        # a1 - Finite element location
        def rule_eq_a1(b, i):
            return blk.ll[i] == b.fe_set[i]*b.Lb
        blk.eq_a1 = Constraint(blk.l, rule=rule_eq_a1)

        # a3 - Reactor Cross-Sectional Area
        blk.eq_a3 = Constraint(expr = blk.Areact == 0.25*blk.pi*blk.Dt**2)
        
        # a4 - Distributor Design
        blk.eq_a4 = Constraint(expr = 10 == blk.nor*blk.Ao*10)

        # ---------------------------------------------------------------------
        # Mass and Energy Balances
        # b1 - Bubble Region Gas Component Balances
        def rule_eq_b1(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dGbcdx[j,i] \
                            - b.dl[i]*b.Ax*b.delta[i]*b.Kbc[j,i] \
                            * (b.cb[j,i]-b.cc[j,i]) + b.Kgbulk_c[j,i]
        blk.eq_b1 = Constraint(blk.l, blk.GasList, rule=rule_eq_b1)

        # b2 - Cloud-Wake Region Gas Component Balances
        def rule_eq_b2(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == b.dl[i]*b.Ax*b.delta[i]*b.Kbc[j,i] \
                        * (b.cb[j,i] - b.cc[j,i]) \
                        - b.dl[i]*b.Ax*b.delta[i]*b.Kce[j,i] \
                        * (b.cc[j,i] - b.ce[j,i]) - b.rgc[j,i]
        blk.eq_b2 = Constraint(blk.l, blk.GasList, rule=rule_eq_b2)

        # b3 - Emulsion Region Gas Component Balances
        def rule_eq_b3(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dGecdx[j,i] \
                        + b.dl[i]*b.Ax*b.delta[i]*b.Kce[j,i] \
                        * (b.cc[j,i] - b.ce[j,i]) - b.Kgbulk_c[j,i] \
                        - b.rge[j,i]
        blk.eq_b3 = Constraint(blk.l, blk.GasList, rule=rule_eq_b3)

        # b4 - Bubble Region Gas Enthalpy Balance
        def rule_eq_b4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dGbhdx[i] \
                            - b.dl[i]*b.Ax*b.delta[i]*b.Hbc[i] \
                            * (b.Tgb[i]-b.Tgc[i]) + b.Hgbulk[i]
        blk.eq_b4 = Constraint(blk.l, rule=rule_eq_b4)

        # b5 - Cloud-Wake Region Gas Enthalpy Balance
        def rule_eq_b5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == b.dl[i]*b.Ax*b.delta[i]*b.Hbc[i] \
                            * (b.Tgb[i]-b.Tgc[i]) \
                        - b.dl[i]*b.Ax*b.delta[i]*b.Hce[i]*(b.Tgc[i]-b.Tge[i])\
                        - b.hp_c[i] \
                        - b.hr_c[i]
        blk.eq_b5 = Constraint(blk.l, rule=rule_eq_b5)

        # b6 - Emulsion Region Gas Enthalpy Balance
        def rule_eq_b6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dGehdx[i] \
                        + b.dl[i]*b.Ax*b.delta[i]*b.Hce[i]*(b.Tgc[i]-b.Tge[i])\
                        - b.hp_e[i] \
                        - b.Hgbulk[i] \
                        - b.hr_e[i]
        blk.eq_b6 = Constraint(blk.l, rule=rule_eq_b6)

        # b7 - Cloud-Wake Region Sorbent Balance
        def rule_eq_b7(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dJcdx[i]*b.Ax + b.Ksbulk[i]
        blk.eq_b7 = Constraint(blk.l, rule=rule_eq_b7)

        # b8 - Emulsion Region Sorbent Balance
        def rule_eq_b8(b, i):
            if (i == 0) or (i == b.nfe):
                return b.Je[i] == b.Jc[i]
            else:
                if b.s_inlet == 'Top' and b.s_outlet == 'Underflow':
                    return b.Je[i]*b.Ax == b.Jc[i]*b.Ax + b.Solid_In_F
                elif b.s_inlet == 'Bottom' and b.s_outlet == 'Overflow':
                    return b.Je[i]*b.Ax == b.Jc[i]*b.Ax - b.Solid_In_F
                else:
                    return b.Je[i] == b.Jc[i]
        blk.eq_b8 = Constraint(blk.l, rule=rule_eq_b8)

        # b9 - Cloud-Wake Region Adsorbed Species Balance
        def rule_eq_b9(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dJccdx[j,i]*b.Ax \
                            + b.Ksbulk_c[j,i] \
                            + b.rsc[j,i] - b.Kcebs_c[j,i]
        blk.eq_b9 = Constraint(blk.l, blk.SolidList, rule=rule_eq_b9)

        # b10 - Emulsion Region Adsorbed Species Balance
        def rule_eq_b10(b, i, j):
            if (i == b.nfe) or (i <= b.feop[0]):
                return Constraint.Skip
            else:
                return 0 == b.dJecdx[j,i]*b.Ax \
                            - b.Ksbulk_c[j,i] \
                            + b.rse[j,i] + b.Kcebs_c[j,i]
        blk.eq_b10 = Constraint(blk.l, blk.SolidList, rule=rule_eq_b10)

        # b11 - Cloud-Wake Region Solid Enthalpy Balance
        def rule_eq_b11(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dJchdx[i]*b.Ax \
                            + b.Hsbulk[i] \
                            - b.Hcebs[i] \
                            + b.hp_c[i] \
                            + b.hr_c[i]
        blk.eq_b11 = Constraint(blk.l, rule=rule_eq_b11)

        # b12 - Emulsion Region Solid Enthalpy Balance
        def rule_eq_b12(b, i):
            if (i == b.nfe) or (i <= b.feop[0]):
                return Constraint.Skip
            else:
                return 0 == b.dJehdx[i]*b.Ax \
                            - b.Hsbulk[i] \
                            + b.Hcebs[i] \
                            + b.hp_e[i] \
                            + b.Qhx[i] \
                            + b.hr_e[i]
        blk.eq_b12 = Constraint(blk.l, rule=rule_eq_b12)

        # ---------------------------------------------------------------------
        # Flowrate and Flux Relationships
        # c1 - Superficial Gas Velocity
        def rule_eq_c1(b, i):
            if i == 0:
                return b.vg[0]*b.Ax == b.Gas_In_F*b.gas_prop_b[0].V_vap
            else:
                return b.vg[i] == b.vb[i]*b.delta[i] + b.ve[i]*(1-b.delta[i])
        blk.eq_c1 = Constraint(blk.l, rule=rule_eq_c1)

        # c2 - Initial emulsion region gas flow rate
        # Assume excess gas goes to bubble region - vmf evaluated at point 1
        blk.eq_c2 = Constraint(expr = blk.sol_prop_e[blk.feop[0]].v_mf*blk.Ax \
                                        == blk.Ge[0]*blk.gas_prop_b[0].V_vap)

        # c3 - Initial bubble region gas flow rate
        blk.eq_c3 = Constraint(expr = blk.Gas_In_F == blk.Ge[0] + blk.Gb[0])

        # c4 - Bubble Gas Velocity
        def rule_eq_c4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Ax*b.delta[i]*b.vb[i] == b.Gb[i]*b.gas_prop_b[i].V_vap
        blk.eq_c4 = Constraint(blk.l, rule=rule_eq_c4)
        
        # c5 - Emulsion Gas Flowrate
        def rule_eq_c5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.ve[i]*b.Ax*b.cet[i]*b.delta_e[i] \
                                == b.Ge[i]
        blk.eq_c5 = Constraint(blk.l, rule=rule_eq_c5)

        # c6 - Bubble Gas Component Flowrate
        def rule_eq_c6(b, i, j):
            return b.Gbc[j,i] == b.Gb[i]*b.yb[j,i]
        blk.eq_c6 = Constraint(blk.l, blk.GasList, rule=rule_eq_c6)

        # c7 - Emulsion Gas Component Flowrate
        def rule_eq_c7(b, i, j):
            return b.Gec[j,i] == b.Ge[i]*b.ye[j,i]
        blk.eq_c7 = Constraint(blk.l, blk.GasList, rule=rule_eq_c7)

        # c8 - Bubble Gas Enthalpy Flowrate
        def rule_eq_c8(b, i):
            return b.Gbh[i] == b.Gb[i]*b.gas_prop_b[i].h_vap
        blk.eq_c8 = Constraint(blk.l, rule=rule_eq_c8)

        # c9 - Emulsion Gas Enthalpy Flowrate
        def rule_eq_c9(b, i):
            if i == 0:
                return b.Geh[i] == b.Ge[i]*b.gas_prop_b[i].h_vap
            else:
                return b.Geh[i] == b.Ge[i]*b.gas_prop_e[i].h_vap
        blk.eq_c9 = Constraint(blk.l, rule=rule_eq_c9)

        # c10 - Cloud-Wake Region Solids Flux
        def rule_eq_c10(b, i):
            if i == 0:
                return b.Jc[i] == b.fw*b.delta[b.feop[i]] \
                                        * b.sol_prop_c[b.feop[i]].rho_sol \
                                        * (1-b.ed[b.feop[i]])*b.vb[b.feop[i]]
            else:
                return b.Jc[i] == b.fw*b.delta[i]*b.sol_prop_c[i].rho_sol \
                                        * (1-b.ed[i])*b.vb[i]
        blk.eq_c10 = Constraint(blk.l, rule=rule_eq_c10)

        # c11 - Emulsion Region Solids Velocity
        def rule_eq_c11(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Je[i] == b.delta_e[i] \
                                * b.sol_prop_e[i].rho_sol*(1-b.ed[i])*b.us[i]
        blk.eq_c11 = Constraint(blk.l, rule=rule_eq_c11)

        # c12 - CW Solids Adsorbed Species Flux
        def rule_eq_c12(b, i, j):
            return b.Jcc[j,i] == b.Jc[i]*b.xc[j,i]
        blk.eq_c12 = Constraint(blk.l, blk.SolidList, rule=rule_eq_c12)

        # c13 - Emulsion Solids Adsorbed Species Flux
        def rule_eq_c13(b, i, j):
            return b.Jec[j,i] == b.Je[i]*b.xe[j,i]
        blk.eq_c13 = Constraint(blk.l, blk.SolidList, rule=rule_eq_c13)

        # c14 - CW Solids Enthalpy Flux
        def rule_eq_c14(b, i):
            if i == 0:
                return b.Tsc[i] == b.Tse[i]
            else:
                return b.Jch[i] == b.Jc[i]*b.sol_prop_c[i].h_sol
        blk.eq_c14 = Constraint(blk.l, rule=rule_eq_c14)

        # c15 - Emulsion Solids Enthalpy Flux
        def rule_eq_c15(b, i):
            if i == b.nfe:
                return b.Tse[i] == b.Tsc[i]
            else:
                return b.Jeh[i] == b.Je[i]*b.sol_prop_e[b.feop[i]].h_sol
        blk.eq_c15 = Constraint(blk.l, rule=rule_eq_c15)

        # ---------------------------------------------------------------------
        # Mole Fraction Relationships
        # e1 - Bubble Region Component Concentrations
        def rule_eq_e1(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.cb[j, i]*b.gas_prop_b[i].V_vap == b.yb[j, i]
        blk.eq_e1 = Constraint(blk.l, blk.GasList, rule=rule_eq_e1)

        # e2 - Cloud-Wake Region Component Concentrations
        def rule_eq_e2(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.cc[j, i] == b.yc[j, i]*b.cct[i]
        blk.eq_e2 = Constraint(blk.l, blk.GasList, rule=rule_eq_e2)

        # e3 - Emulsion Region Component Concentrations
        def rule_eq_e3(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.ce[j, i] == b.ye[j, i]*b.cet[i]
        blk.eq_e3 = Constraint(blk.l, blk.GasList, rule=rule_eq_e3)

        # e4 - Bubble Region Total Concentration
        def rule_eq_e4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 1 == sum(b.cb[j, i] for j in b.GasList)\
                                * b.gas_prop_b[i].V_vap
        blk.eq_e4 = Constraint(blk.l, rule=rule_eq_e4)

        # e5 - Cloud-Wake Region Total Concentration
        def rule_eq_e5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.cct[i] == sum(b.cc[j, i] for j in b.GasList)
        blk.eq_e5 = Constraint(blk.l, rule=rule_eq_e5)

        # e6 - Emulsion Region Total Concentration
        def rule_eq_e6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return blk.cet[i] == sum(blk.ce[j, i] for j in blk.GasList)
        blk.eq_e6 = Constraint(blk.l, rule=rule_eq_e6)

        # ---------------------------------------------------------------------
        # Bulk Flow and Mixing Relationships
        # f1 - Bulk Gas Mass Transfer
        def rule_eq_f1(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Kgbulk_c[j,i]*b.db[i] == (6*b.Kd*b.delta[i] \
                                * b.dl[i]*b.Ax)\
                                * (0.5*((b.cet[i]-(1/b.gas_prop_b[i].V_vap)) \
                                + sqrt((b.cet[i]-(1/b.gas_prop_b[i].V_vap))**2\
                                    + b.eps**2))*b.yb[j,i] \
                                - 0.5*(-(b.cet[i]-(1/b.gas_prop_b[i].V_vap)) \
                                + sqrt((b.cet[i]-(1/b.gas_prop_b[i].V_vap))**2\
                                    + b.eps**2))*b.ye[j,i])
        blk.eq_f1 = Constraint(blk.l, blk.GasList, rule=rule_eq_f1)

        # f2 - Bulk Gas Enthalpy Transfer
        def rule_eq_f2(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Hgbulk[i]*b.db[i] == (6*b.Kd*b.delta[i]*b.dl[i]*b.Ax)\
                                * (0.5*((b.cet[i]-(1/b.gas_prop_b[i].V_vap)) \
                                + sqrt((b.cet[i]-(1/b.gas_prop_b[i].V_vap))**2\
                                    + b.eps**2))*b.gas_prop_e[i].h_vap \
                                - 0.5*(-(b.cet[i]-(1/b.gas_prop_b[i].V_vap)) \
                                + sqrt((b.cet[i]-(1/b.gas_prop_b[i].V_vap))**2\
                                    + b.eps**2))*b.gas_prop_b[i].h_vap)
        blk.eq_f2 = Constraint(blk.l, rule=rule_eq_f2)

        # f3 - Bulk Solids Mass Transfer
        def rule_eq_f3(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Ksbulk_c[j,i] == 0.5*(b.Ksbulk[i] \
                                        + sqrt(b.Ksbulk[i]**2 + b.eps**2)) \
                                            * b.xe[j,b.feom[i]] \
                                    - 0.5*(-b.Ksbulk[i] \
                                        + sqrt(b.Ksbulk[i]**2 + b.eps**2)) \
                                            * b.xc[j,i]
        blk.eq_f3 = Constraint(blk.l, blk.SolidList, rule=rule_eq_f3)

        # f4 - Bulk Solids Enthalpy Transfer
        def rule_eq_f4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Hsbulk[i] == 0.5*(b.Ksbulk[i] \
                                        + sqrt(b.Ksbulk[i]**2 + b.eps**2)) \
                                            * b.sol_prop_e[i].h_sol \
                                    - 0.5*(-b.Ksbulk[i] \
                                        + sqrt(b.Ksbulk[i]**2 + b.eps**2)) \
                                            * b.sol_prop_c[i].h_sol
        blk.eq_f4 = Constraint(blk.l, rule=rule_eq_f4)

        # f5 - Bulk Solids Component Mixing
        def rule_eq_f5(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Kcebs_c[j,i] == b.dl[i]*b.Ax*b.delta[i] \
                                            * b.sol_prop_e[i].rho_sol \
                                            * b.Kcebs[i]*(b.xc[j,i] \
                                                - b.xe[j,b.feom[i]])
        blk.eq_f5 = Constraint(blk.l, blk.SolidList, rule=rule_eq_f5)
        
        # f6 - Bulk Solids Enthalpy Mixing
        def rule_eq_f6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Hcebs[i] == b.dl[i]*b.Ax*b.delta[i]\
                                            * b.sol_prop_e[i].rho_sol \
                                            * b.Kcebs[i] \
                                            * (b.sol_prop_c[i].h_sol \
                                               - b.sol_prop_e[i].h_sol)
        blk.eq_f6 = Constraint(blk.l, rule=rule_eq_f6)
        
        # f7 - Cloud-Wake Region Gas-Solids Convective Enthalpy Transfer
        def rule_eq_f7(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hp_c[i] == b.dl[i]*b.Ax*b.delta[i]*b.fcw[i] \
                            * (1-b.ed[i]) \
                            * b.sol_prop_c[i].rho_sol \
                            * b.sol_prop_c[i].ap*b.hp[i]*(b.Tgc[i]-b.Tsc[i])
        blk.eq_f7 = Constraint(blk.l, rule=rule_eq_f7)
        
        # f8 - Emulsion Region Gas-Solids Convective Enthalpy Transfer
        def rule_eq_f8(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hp_e[i] == b.dl[i]*b.Ax \
                            * b.delta_e[i] \
                            * (1-b.ed[i])*b.sol_prop_e[i].rho_sol \
                            * b.sol_prop_e[i].ap*b.hp[i] \
                            * (b.Tge[i]-b.Tse[b.feom[i]])
        blk.eq_f8 = Constraint(blk.l, rule=rule_eq_f8)
        
        # f9 - Cloud-Wake Region Gas-Solids Enthalpy Transfer from Reaction
        def rule_eq_f9(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hr_c[i] == sum(b.rgc[j,i] for j in b.GasList) \
                            * b.gas_prop_c[i].h_vap
        blk.eq_f9 = Constraint(blk.l, rule=rule_eq_f9)
        
        # f10 - Emulsion Region Gas-Solids Enthalpy Transfer from Reaction
        def rule_eq_f10(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hr_e[i] == sum(b.rge[j,i] for j in b.GasList) \
                            * b.gas_prop_e[i].h_vap
        blk.eq_f10 = Constraint(blk.l, rule=rule_eq_f10)

        # ---------------------------------------------------------------------
        # g1 - Pressure Drop
        def rule_eq_g1(b, i):
            if i == 0:
                return b.P[i] == b.Gas_In_P - 3400
            else:
                return b.dPdx[i] == - b.dl[i]*(1-b.e[i]) \
                                            * b.sol_prop_e[i].rho_sol*b.gc
        blk.eq_g1 = Constraint(blk.l, rule=rule_eq_g1)

        # h1 - Archimedes Number
        def rule_eq_h1(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Ar[i]*b.gas_prop_e[i].mu_vap**2 \
                            == (b.sol_prop_e[i].dp**3)*b.gas_prop_e[i].rho_vap\
                                * (b.sol_prop_e[i].rho_sol \
                                    -b.gas_prop_e[i].rho_vap)*b.gc
        blk.eq_h1 = Constraint(blk.l, rule=rule_eq_h1)
        
        # h12 - Emulsion Region Volume Fraction
        def rule_eq_h12(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.delta_e[i] == 1 - (b.fcw[i] + 1)*b.delta[i]
        blk.eq_h12 = Constraint(blk.l, rule=rule_eq_h12)

    def _make_bdry_conds(self):
        """
        Create solids boundary conditions (emulsion region).
        """
        blk = self
        
        # d1-d4 - Gas inlet condiitons
        def rule_d1(b, j):
            return b.yb[j,0] == b.Gas_In_y[j]
        blk.eq_d1 = Constraint(blk.GasList, rule=rule_d1)
        def rule_d2(b, j):
            return b.ye[j,0] == b.Gas_In_y[j]
        blk.eq_d2 = Constraint(blk.GasList, rule=rule_d2)
        blk.eq_d3 = Constraint(expr = blk.Tgb[0] == blk.Gas_In_T)
        blk.eq_d4 = Constraint(expr = blk.Tge[0] == blk.Gas_In_T)

        # ---------------------------------------------------------------------
        # d5-d8 - Gas outlet conditions
        blk.eq_d5 = Constraint(expr = blk.Gas_Out_F == blk.Gb[blk.nfe] \
                                                        + blk.Ge[blk.nfe])
        blk.eq_d6 = Constraint(expr = blk.Gas_Out_F*blk.gas_prop_out.h_vap \
                            == blk.Gb[blk.nfe]*blk.gas_prop_b[blk.nfe].h_vap \
                            + blk.Ge[blk.nfe]*blk.gas_prop_e[blk.nfe].h_vap)
        blk.eq_d7 = Constraint(expr = blk.Gas_Out_P == blk.P[blk.nfe])
        def rule_d8(b, j):
            return b.Gas_Out_F*b.Gas_Out_y[j] == b.Gb[b.nfe]*b.yb[j,b.nfe] \
                                                    + b.Ge[b.nfe]*b.ye[j,b.nfe]
        blk.eq_d8 = Constraint(blk.GasList, rule=rule_d8)
        
        # ---------------------------------------------------------------------
        # Solid Recirculation Boundary Conditions
        def rule_eq_b17(b, j):
            return b.Jcc[j,0] == b.Jec[j,0]
        blk.eq_b17 = Constraint(blk.SolidList, rule=rule_eq_b17)
        
        def rule_eq_b18(b, j):
            return b.Jec[j,b.nfe] == b.Jcc[j,b.nfe]
        blk.eq_b18 = Constraint(blk.SolidList, rule=rule_eq_b18)
        
        blk.eq_b19 = Constraint(expr = blk.Jch[0] == blk.Jeh[0])
        blk.eq_b20 = Constraint(expr = blk.Jeh[blk.nfe] == blk.Jch[blk.nfe])

        # ---------------------------------------------------------------------
        if self.s_inlet == 'Bottom' and self.s_outlet == 'Underflow':
        # Bottom Feed, Underflow Outlet
            # Adsorbed Species Balances
            def rule_eq_b13(b, j):
                return 0 == b.dJecdx[j,b.feop[0]]*b.Ax \
                        + b.Solid_In_F*b.Solid_In_x[j] \
                        - b.Solid_Out_F*b.Solid_Out_x[j] \
                        - b.Ksbulk_c[j,b.feop[0]] \
                        + b.rse[j,b.feop[0]] + b.Kcebs_c[j,b.feop[0]]
            blk.eq_b13 = Constraint(blk.SolidList, rule=rule_eq_b13)

            def rule_eq_b14(b, j):
                return 0 == b.dJecdx[j,b.nfe]*b.Ax \
                        - b.Ksbulk_c[j,b.nfe] + b.rse[j,b.nfe] \
                        + b.Kcebs_c[j,b.nfe]
            blk.eq_b14 = Constraint(blk.SolidList, rule=rule_eq_b14)

            # Enthalpy Balances
            blk.eq_b15 = Constraint(expr = 0 == blk.dJehdx[blk.feop[0]]*blk.Ax\
                    + blk.Solid_In_F*blk.sol_prop_f.h_sol \
                    - blk.Solid_Out_F*blk.sol_prop_e[blk.feop[0]].h_sol \
                    - blk.Hsbulk[blk.feop[0]] \
                    + blk.Hcebs[blk.feop[0]]\
                    + blk.hp_e[blk.feop[0]] \
                    + blk.Qhx[blk.feop[0]] \
                    + blk.hr_e[blk.feop[0]])

            blk.eq_b16 = Constraint(expr = 0 == blk.dJehdx[blk.nfe]*blk.Ax \
                    - blk.Hsbulk[blk.nfe] \
                    + blk.Hcebs[blk.nfe] \
                    + blk.hp_e[blk.nfe] \
                    + blk.Qhx[blk.nfe] \
                    + blk.hr_e[blk.nfe])

            # d9-d11 - Solid outlets
            blk.eq_d9 = Constraint(expr = blk.Solid_Out_F == blk.Solid_In_F)
            blk.eq_d10 = Constraint(expr = blk.Solid_Out_T \
                                            == blk.Tse[0])
            def rule_eq_d11(b, j):
                return b.Solid_Out_x[j] == b.xe[j,0]
            blk.eq_d11 = Constraint(blk.SolidList, rule=rule_eq_d11)

        # ---------------------------------------------------------------------
        elif self.s_inlet == 'Top' and self.s_outlet == 'Underflow':
        # Top Feed, Underflow Outlet
            # Adsorbed Species Balances
            def rule_eq_b13(b, j):
                return 0 == b.dJecdx[j,b.feop[0]]*b.Ax \
                        - b.Solid_Out_F*b.Solid_Out_x[j] \
                        - b.Ksbulk_c[j,b.feop[0]] \
                        + b.rse[j,b.feop[0]] + b.Kcebs_c[j,b.feop[0]]
            blk.eq_b13 = Constraint(blk.SolidList, rule=rule_eq_b13)

            def rule_eq_b14(b, j):
                return 0 == b.dJecdx[j,b.nfe]*b.Ax \
                        + b.Solid_In_F*b.Solid_In_x[j] \
                        - b.Ksbulk_c[j,b.nfe] + b.rse[j,b.nfe] \
                        + b.Kcebs_c[j,b.nfe]
            blk.eq_b14 = Constraint(blk.SolidList, rule=rule_eq_b14)

            # Enthalpy Balances
            blk.eq_b15 = Constraint(expr = 0 == blk.dJehdx[blk.feop[0]]*blk.Ax\
                    - blk.Solid_Out_F*blk.sol_prop_e[blk.feop[0]].h_sol \
                    - blk.Hsbulk[blk.feop[0]] \
                    + blk.Hcebs[blk.feop[0]]\
                    + blk.hp_e[blk.feop[0]] \
                    + blk.Qhx[blk.feop[0]] \
                    + blk.hr_e[blk.feop[0]])

            blk.eq_b16 = Constraint(expr = 0 == blk.dJehdx[blk.nfe]*blk.Ax \
                    + blk.Solid_In_F*blk.sol_prop_f.h_sol \
                    - blk.Hsbulk[blk.nfe] \
                    + blk.Hcebs[blk.nfe] \
                    + blk.hp_e[blk.nfe] \
                    + blk.Qhx[blk.nfe] \
                    + blk.hr_e[blk.nfe])

            # d9-d11 - Solid outlets
            blk.eq_d9 = Constraint(expr = blk.Solid_Out_F == blk.Solid_In_F)
            blk.eq_d10 = Constraint(expr = blk.Solid_Out_T \
                                            == blk.Tse[0])
            def rule_eq_d11(b, j):
                return b.Solid_Out_x[j] == b.xe[j,0]
            blk.eq_d11 = Constraint(blk.SolidList, rule=rule_eq_d11)

        # ---------------------------------------------------------------------
        elif self.s_inlet == 'Bottom' and self.s_outlet == 'Overflow':
        # Bottom Feed, Overflow Outlet
            # Adsorbed Species Balances
            def rule_eq_b13(b, j):
                return 0 == b.dJecdx[j,b.feop[0]]*b.Ax \
                        + b.Solid_In_F*b.Solid_In_x[j] \
                        - b.Ksbulk_c[j,b.feop[0]] \
                        + b.rse[j,b.feop[0]] + b.Kcebs_c[j,b.feop[0]]
            blk.eq_b13 = Constraint(blk.SolidList, rule=rule_eq_b13)

            def rule_eq_b14(b, j):
                return 0 == b.dJecdx[j,b.nfe]*b.Ax \
                        - b.Solid_Out_F*b.Solid_Out_x[j] \
                        - b.Ksbulk_c[j,b.nfe] + b.rse[j,b.nfe] \
                        + b.Kcebs_c[j,b.nfe]
            blk.eq_b14 = Constraint(blk.SolidList, rule=rule_eq_b14)

            # Enthalpy Balances
            blk.eq_b15 = Constraint(expr = 0 == blk.dJehdx[blk.feop[0]]*blk.Ax\
                    + blk.Solid_In_F*blk.sol_prop_f.h_sol \
                    - blk.Hsbulk[blk.feop[0]] \
                    + blk.Hcebs[blk.feop[0]]\
                    + blk.hp_e[blk.feop[0]] \
                    + blk.Qhx[blk.feop[0]] \
                    + blk.hr_e[blk.feop[0]])

            blk.eq_b16 = Constraint(expr = 0 == blk.dJehdx[blk.nfe]*blk.Ax \
                    - blk.Solid_Out_F*blk.sol_prop_e[blk.feom[blk.nfe]].h_sol \
                    - blk.Hsbulk[blk.nfe] \
                    + blk.Hcebs[blk.nfe] \
                    + blk.hp_e[blk.nfe] \
                    + blk.Qhx[blk.nfe] \
                    + blk.hr_e[blk.nfe])

            # d9-d11 - Solid outlets
            blk.eq_d9 = Constraint(expr = blk.Solid_Out_F == blk.Solid_In_F)
            blk.eq_d10 = Constraint(expr = blk.Solid_Out_T \
                                            == blk.Tse[blk.feom[blk.nfe]])
            def rule_eq_d11(b, j):
                return b.Solid_Out_x[j] == b.xe[j,blk.feom[blk.nfe]]
            blk.eq_d11 = Constraint(blk.SolidList, rule=rule_eq_d11)

        # ---------------------------------------------------------------------
        else:
        # Top Feed, Overflow Outlet
            # Adsorbed Species Balances
            def rule_eq_b13(b, j):
                return 0 == b.dJecdx[j,b.feop[0]]*b.Ax \
                        - b.Ksbulk_c[j,b.feop[0]] \
                        + b.rse[j,b.feop[0]] + b.Kcebs_c[j,b.feop[0]]
            blk.eq_b13 = Constraint(blk.SolidList, rule=rule_eq_b13)

            def rule_eq_b14(b, j):
                return 0 == b.dJecdx[j,b.nfe]*b.Ax \
                        + b.Solid_In_F*b.Solid_In_x[j] \
                        - b.Solid_Out_F*b.Solid_Out_x[j] \
                        - b.Ksbulk_c[j,b.nfe] + b.rse[j,b.nfe] \
                        + b.Kcebs_c[j,b.nfe]
            blk.eq_b14 = Constraint(blk.SolidList, rule=rule_eq_b14)

            # Enthalpy Balances
            blk.eq_b15 = Constraint(expr = 0 == blk.dJehdx[blk.feop[0]]*blk.Ax\
                    - blk.Hsbulk[blk.feop[0]] \
                    + blk.Hcebs[blk.feop[0]]\
                    + blk.hp_e[blk.feop[0]] \
                    + blk.Qhx[blk.feop[0]] \
                    + blk.hr_e[blk.feop[0]])

            blk.eq_b16 = Constraint(expr = 0 == blk.dJehdx[blk.nfe]*blk.Ax \
                    + blk.Solid_In_F*blk.sol_prop_f.h_sol \
                    - blk.Solid_Out_F*blk.sol_prop_e[blk.feom[blk.nfe]].h_sol \
                    - blk.Hsbulk[blk.nfe] \
                    + blk.Hcebs[blk.nfe] \
                    + blk.hp_e[blk.nfe] \
                    + blk.Qhx[blk.nfe] \
                    + blk.hr_e[blk.nfe])

            # d9-d11 - Solid outlets
            blk.eq_d9 = Constraint(expr = blk.Solid_Out_F == blk.Solid_In_F)
            blk.eq_d10 = Constraint(expr = blk.Solid_Out_T \
                                            == blk.Tse[blk.feom[blk.nfe]])
            def rule_eq_d11(b, j):
                return b.Solid_Out_x[j] == b.xe[j,blk.feom[blk.nfe]]
            blk.eq_d11 = Constraint(blk.SolidList, rule=rule_eq_d11)

    def _make_hydro_vars(self):
        """
        Create the variables needed for the hydrodynamic model.
        
        Group: H, I
        """
        blk = self

        # Create variables internal to hydrodynamic model
        blk.dbm = Var(blk.l, domain=Reals,
                    doc='Maximum Theoretical Bubble Diameter [m]')
        blk.g1 = Var(blk.l, domain=Reals,
                    doc='Bubble Growth Coefficient [-]')
        blk.vbr = Var(blk.l, domain=Reals,
                    doc='Bubble Rise Velocity [m/s]')
        
        blk.ddbdx = DerivativeVar(blk.db, wrt=blk.l)

    def _make_hydro_model(self):
        """
        Create the hydrodynamic model for the bed.
        
        Group: H, I
        """     
        blk = self

        # Fluidisation Conditions and Bubble Behaviour
        # h2 - Emulsion Region Gas Velocity
        def rule_eq_h2(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.ve[i]*((b.sol_prop_e[i].dp**0.568)*(b.gc**0.663) \
                            * ((b.sol_prop_e[i].rho_sol \
                                - b.gas_prop_e[i].rho_vap)**0.663) \
                            * ((b.ll[i])**0.244)) \
                        == b.sol_prop_e[i].v_mf*188 \
                            * (b.gas_prop_e[i].rho_vap**0.089) \
                            * (b.gas_prop_e[i].mu_vap**0.371) \
                            * exp(0.508*b.sol_prop_e[i].F)
        blk.eq_h2 = Constraint(blk.l, rule=rule_eq_h2)

        # h3 - Average Cross-Sectional Voidage
        def rule_eq_h3(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return (1-b.e[i]) == (1-b.delta[i])*(1-b.ed[i])
        blk.eq_h3 = Constraint(blk.l, rule=rule_eq_h3)

        # h4 - Bubble Size Coefficients
        def rule_eq_h4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return (blk.g1[i]*blk.sol_prop_e[i].v_mf)**2 \
                                        == (2.56E-2**2)*(blk.Dt/blk.gc)
        blk.eq_h4 = Constraint(blk.l, rule=rule_eq_h4)

        # h5 - Maximum Bubble Diameter
        def rule_eq_h5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return (b.dbm[i]**5)*b.gc \
                                    == (2.59**5)*((b.vg[i]-b.ve[i])*b.Ax)**2
        blk.eq_h5 = Constraint(blk.l, rule=rule_eq_h5)

        # h6 - Constrained Bubble Diameter
        def rule_eq_h6(b, i):
            if i == 0:
                return b.db[0] == 1.38*(b.gc**(-0.2)) \
                                * ((b.vg[0]-b.ve[b.feop[0]])*b.Ao)**0.4
            else:
                return b.ddbdx[i]*b.Dt == 0.3*b.dl[i]*(b.dbm[i] \
                                        - b.db[i] - b.g1[i] \
                                            * (b.Dt*b.db[i])**0.5)
        blk.eq_h6 = Constraint(blk.l, rule=rule_eq_h6)

        # h7 - Bubble Rise Velocity
        def rule_eq_h7(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.vbr[i]**2 == (0.711**2)*(b.gc*b.db[i])
        blk.eq_h7 = Constraint(blk.l, rule=rule_eq_h7)

        # h8 - Bubble Velocity
        def rule_eq_h8(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                if b.vb_method == 'Werther A':
                    return b.vb[i] == 1.55*((b.vg[i]-b.sol_prop_e[i].v_mf) \
                            + 14.1*(b.db[i]+0.005))*(b.Dte**0.32) + b.vbr[i]
                elif b.vb_method == 'Werther B':
                    return b.vb[i] == 1.6*((b.vg[i]-b.sol_prop_e[i].v_mf) \
                            + 1.13*(b.db[i]**0.5))*(b.Dte**1.35) + b.vbr[i]
                else:               # Davidson Model
                    return b.vb[i] == b.vg[i] - b.sol_prop_e[i].v_mf + b.vbr[i]
        blk.eq_h8 = Constraint(blk.l, rule=rule_eq_h8)

        # h9 - Cloud to Bubble Volume Ratio
        def rule_eq_h9(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 3*b.sol_prop_e[i].v_mf/b.sol_prop_e[i].e_mf \
                        == b.fc[i]*(b.vbr[i]-b.sol_prop_e[i].v_mf \
                            / b.sol_prop_e[i].e_mf)
                '''return 6*(b.sol_prop_e[i].v_mf/b.sol_prop_e[i].e_mf) \
                        == b.fc[i]*(b.vbr[i]-(b.sol_prop_e[i].v_mf \
                            / b.sol_prop_e[i].e_mf)*(1 - b.sol_prop_e[i].v_mf\
                                / b.sol_prop_e[i].e_mf) \
                        + (((1 + b.sol_prop_e[i].v_mf/b.sol_prop_e[i].e_mf) \
                            * b.sol_prop_e[i].v_mf/b.sol_prop_e[i].e_mf \
                                - b.vbr[i])**2 \
                        + b.eps**2)**0.5)'''
        blk.eq_h9 = Constraint(blk.l, rule=rule_eq_h9)

        # h10 - Cloud-Wake to Bubble Volume Ratio
        def rule_eq_h10(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.fcw[i] == b.fc[i] + b.fw
        blk.eq_h10 = Constraint(blk.l, rule=rule_eq_h10)

        # h11 - Emulsion Region Voidage
        def rule_eq_h11(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return (1 - b.sol_prop_e[i].e_mf)*((b.sol_prop_e[i].dp**0.1) \
                        * (b.gc**0.118)*((b.sol_prop_e[i].rho_sol \
                            - b.gas_prop_e[i].rho_vap)**0.118) \
                            * ((b.ll[i])**0.043)) \
                        == 2.54*(b.gas_prop_e[i].rho_vap**0.016) \
                            * (b.gas_prop_e[i].mu_vap**0.066) \
                            * exp(0.090*b.sol_prop_e[i].F)*(1-b.ed[i])
        blk.eq_h11 = Constraint(blk.l, rule=rule_eq_h11)

        # ---------------------------------------------------------------------
        # K-L Model Heat and Mass Transfer Coefficients
        # i1 - Bubble to Cloud-Wake Gas Mass Transfer Coefficient
        def rule_eq_i1(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Kbc[j,i]*b.db[i]**(5/4) \
                    == 1.32*4.5*b.sol_prop_c[i].v_mf*b.db[i]**(1/4) \
                        + 5.85*b.gas_prop_e[i].D_vap[j]**0.5*b.gc**(1/4)
        blk.eq_i1 = Constraint(blk.l, blk.GasList, rule=rule_eq_i1)

        # i2 - Cloud-Wake to Emuslion Gas Mass Transfer Coefficient
        def rule_eq_i2(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Kce[j,i]**2*(b.db[i]**3) == (6.77**2)*b.ed[i] \
                                            * b.gas_prop_e[i].D_vap[j]\
                                            * b.vbr[i]
        blk.eq_i2 = Constraint(blk.l, blk.GasList, rule=rule_eq_i2)

        # i3 - Cloud-Wake to Emulsion Solid Mass Transfer Coefficient
        def rule_eq_i3(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Kcebs[i]*((1-b.delta[i])*b.ed[i]*b.db[i]) \
                            == 3*(1-b.ed[i])*b.ve[i]
        blk.eq_i3 = Constraint(blk.l, rule=rule_eq_i3)

        # i4 - Bubble to Cloud-Wake Gas Heat Transfer Coefficient
        def rule_eq_i4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Hbc[i]*b.db[i]**(5/4) == 1.32*4.5 \
                                    * b.sol_prop_c[i].v_mf\
                                    * b.gas_prop_b[i].cp_vap*b.db[i]**(1/4) \
                                    / b.gas_prop_b[i].V_vap \
                                + 5.85*(b.gas_prop_b[i].k_vap \
                                    * b.gas_prop_b[i].cp_vap \
                                    / b.gas_prop_b[i].V_vap)**0.5*(b.gc**0.25)
        blk.eq_i4 = Constraint(blk.l, rule=rule_eq_i4)

        # i5 - Cloud-Wake to Emulsion Gas Heat Transfer Coefficient
        def rule_eq_i5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Hce[i]**2*b.db[i]**3 == 6.78**2*b.ed[i]*b.vbr[i] \
                                    * b.gas_prop_e[i].k_vap*b.cct[i] \
                                    * b.gas_prop_e[i].cp_vap
        blk.eq_i5 = Constraint(blk.l, rule=rule_eq_i5)

        # i6 - Convective Heat Transfer Coefficient
        def rule_eq_i6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hp[i]*b.sol_prop_e[i].dp == 0.03 \
                                    * b.gas_prop_e[i].k_vap \
                                    * (b.ve[i]*b.sol_prop_e[i].dp \
                                    * b.gas_prop_e[i].rho_vap \
                                        / b.gas_prop_e[i].mu_vap)**1.3
        blk.eq_i6 = Constraint(blk.l, rule=rule_eq_i6)

    def _make_HX_vars(self):
        """
        Create the variables needed for the heat exchanger.
        
        Group: J, K
        """
        blk = self

        # Create HX tube model
        # Heat Exchanger Dimensions
        blk.Ahx = Var(domain=Reals,
                    doc='Total Area of Heat Exchanger Surfaces [m^2]')
        blk.dx = Var(domain=Reals,
                    doc='Diameter of Heat Exchanger Tubes [m]')                    
        blk.lhx = Var(domain=Reals,
                    doc='Heat Exchanger Tube Spacing (Pitch-Diameter) [m]')
        blk.lp = Var(domain=Reals,
                    doc='Heat Exchanger Tube Pitch [m]')
        blk.Nx = Var(domain=Reals,
                    doc='Number of Heat Exchanger Tubes [-]')
        
        # HX Fluid Inlet Conditions
        blk.HX_In_F = Var(domain=Reals,
                    doc='Heat Exchanger Fluid Flowrate at Inlet [mol/s]')
        blk.HX_In_T = Var(domain=Reals,
                    doc='Heat Exchanger Fluid Temperature at Inlet [K]')
        blk.HX_In_P = Var(domain=Reals,
                    doc='Heat Exchanger Fluid Pressure at Inlet [Pa]')
        blk.HX_In_y = Var(blk.HXList, domain=Reals,
                    doc='Heat Exchanger Fluid Mole Fractions [mol/mol]')

        # HX Fluid State
        blk.Hhx = Var(blk.l, domain=Reals,
                    doc='Heat Exchanger Fluid Enthalpy Flow [J/s]')
        blk.Phx = Var(blk.l, domain=Reals,
                    doc='Heat Exchanger Fluid Pressure [Pa]')
        blk.Thx = Var(blk.l, domain=Reals,
                doc='Heat Exchanger Fluid Temperature [K]')
        blk.Ttube = Var(blk.l, domain=Reals,
                    doc='Heat Exchanger Tube Temperature [K]')
        blk.dThx = Var(blk.l,
                    doc='Temperature Difference between HX Tubes and Bed [K]')

        # HX Tube Heat Transfer Coefficients
        blk.fb = Var(blk.l, domain=Reals,
                    doc='Fraction of Time HX Tubes Contact Dense Packets [-]')
        blk.fn = Var(blk.l, domain=Reals,
                    doc='Fluidisation Number [-]') 
        blk.hd = Var(blk.l, domain=Reals,
                    doc='Convective Heat Transfer Coefficient of Dense '\
                    'Packets [J/m^2.K.s]')
        blk.hl = Var(blk.l, domain=Reals,
                    doc='Convective Heat Transfer Coefficient of Gas Bubbles '\
                    '[J/m^2.K.s]') 
        blk.ht = Var(blk.l, domain=Reals,
                    doc='Overall Convective Heat Transfer Coefficient '\
                    '[J/m^2.K.s]') 
        blk.kpa = Var(blk.l, domain=Reals, initialize=0.1,
                    doc='Thermal Conductivity of Bed at Minimum Fluidisation '\
                    '[J/b.K.s]') 
        blk.Pr = Var(blk.l, domain=Reals, initialize=0.7,
                    doc='Prandlt Number [-]')
        blk.tau = Var(blk.l, domain=Reals, initialize=0.1,
                    doc='Average Residence Time of Dense Packets at '\
                    'HX Surface [s]')
        
        # HX Fluid Inlet Conditions
        blk.HX_Out_F = Var(domain=Reals,
                    doc='HX Fluid Flowrate at Outlet [mol/s]')
        blk.HX_Out_T = Var(domain=Reals,
                    doc='HX Fluid Temperature at Outlet [K]')
        blk.HX_Out_P = Var(domain=Reals,
                    doc='HX Fluid Pressure at Outlet [Pa]')
        blk.HX_Out_y = Var(blk.HXList, domain=Reals,
                    doc='HX Fluid Mole fractions at Outlet [mol/mol]')
        
        # Derivative Variables
        blk.dHhxdx = DerivativeVar(blk.Hhx, wrt=blk.l)
        blk.dPhxdx = DerivativeVar(blk.Phx, wrt=blk.l)

    def _make_HX_model(self):
        """
        Create the heat exchanger model for the bed.
        
        Group: J, K
        """
        blk = self

        # Create HX tube model
        # a5 - Total Reactor Cross-Sectional Area (incl. HX tubes)
        blk.eq_a5 = Constraint(expr = blk.Areact == blk.Ax \
                                                + (blk.pi/4)*blk.dx**2*blk.Nx)

        # a6 - Hydraulic Diameter
        blk.eq_a6 = Constraint(expr = blk.Ax == 0.25*blk.Dte*blk.pi \
                                                * (blk.Dt+blk.dx*blk.Nx))
            
        # a7, a8 - HX Tube Pitch and Spacing
        blk.eq_a7 = Constraint(expr = blk.Areact == blk.Nx*blk.lp**2)
        blk.eq_a8 = Constraint(expr = blk.lhx == blk.lp - blk.dx)

        # a9 - Surface Area of HX Tubes
        blk.eq_a9 = Constraint(expr = blk.Ahx \
                                                == blk.pi*blk.dx*blk.Lb*blk.Nx)

        # -----------------------------------------------------------------
        # Mickley and Fairbanks HX Model
        # j1 - Thermal Conductivity of Bed at Minimum Fluidisation
        def rule_eq_j1(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.kpa[i] == (3.58-2.5*b.ed[i])*b.gas_prop_e[i].k_vap \
                                * ((b.sol_prop_e[i].k_sol \
                                    / b.gas_prop_e[i].k_vap) \
                                ** (0.46-0.46*b.ed[i]))
        blk.eq_j1 = Constraint(blk.l, rule=rule_eq_j1)

        # j2 - Fluidisation Number
        def rule_eq_j2(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.fn[i]*b.sol_prop_e[i].v_mf == b.vg[i]
        blk.eq_j2 = Constraint(blk.l, rule=rule_eq_j2)

        # j3 - Residence Time of Emulsion Packets at HX Surface
        def rule_eq_j3(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.tau[i] == 0.44*((b.sol_prop_e[i].dp*b.gc \
                                / ((b.sol_prop_e[i].v_mf**2) \
                                * ((b.fn[i]-b.ah)**2)))**0.14) \
                                * ((b.sol_prop_e[i].dp/b.dx)**0.225)
        blk.eq_j3 = Constraint(blk.l, rule=rule_eq_j3)

        # j4 - Fraction of Time HX Surface is Exposed to Emulsion Packets
        def rule_eq_j4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.fb[i] == 0.33*(((b.sol_prop_e[i].v_mf**2) \
                                * ((b.fn[i]-b.ah)**2) \
                                / (b.sol_prop_e[i].dp*b.gc))**0.14)
        blk.eq_j4 = Constraint(blk.l, rule=rule_eq_j4)

        # j5 - Dense Region Heat Transfer Coefficient
        def rule_eq_j5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hd[i]*sqrt(b.pi*b.tau[i]) == 2*sqrt(b.kpa[i] \
                                * b.sol_prop_e[i].rho_sol \
                                * b.sol_prop_e[i].cp_sol*(1-b.ed[i]))
        blk.eq_j5 = Constraint(blk.l, rule=rule_eq_j5)

        # j6 - Bubble Region Heat Transfer Coefficient - Prandlt Number
        def rule_eq_j6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Pr[i]*b.gas_prop_e[i].k_vap \
                                == b.gas_prop_e[i].cp_vap \
                                    * b.gas_prop_e[i].mu_vap\
                                    / b.gas_prop_e[i].MW_vap
        blk.eq_j6 = Constraint(blk.l, rule=rule_eq_j6)

        # j7 - Bubble Region Heat Transfer Coefficient
        def rule_eq_j7(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.hl[i]*b.sol_prop_e[i].dp == 0.009*(b.Ar[i]**0.5) \
                                * (b.Pr[i]**0.33)*b.gas_prop_e[i].k_vap
        blk.eq_j7 = Constraint(blk.l, rule=rule_eq_j7)

        # j8 - Total HX Heat Transfer Coefficient
        def rule_eq_j8(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.ht[i] == b.fb[i]*b.hd[i] + (1-b.fb[i])*b.hl[i]
        blk.eq_j8 = Constraint(blk.l, rule=rule_eq_j8)

        # ---------------------------------------------------------------------
        # k1 - Total Heat Duty
        blk.eq_k1 = Constraint(expr = blk.Q == blk.HX_In_F \
                                    * (blk.hx_prop_in.h_mix \
                                       - blk.hx_prop[0].h_mix))

        # k2 - HX Tube Heat Transfer
        def rule_eq_k2(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.Qhx[i] == b.dl[i]*b.pi*b.dx*b.ht[i]*b.dThx[i] \
                                            * b.Nx*b.Cr
        blk.eq_k2 = Constraint(blk.l, rule=rule_eq_k2)

        # k3 - HX Fluid Pressure Drop
        def rule_eq_k3(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.dPhxdx[i] == -b.hx_prop[i].rho_mix*b.gc*b.dl[i]
        blk.eq_k3 = Constraint(blk.l, rule=rule_eq_k3)

        # k4 - HX Fluid Energy Balance
        def rule_eq_k4(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return 0 == -b.dHhxdx[i] + b.Qhx[i]
        blk.eq_k4 = Constraint(blk.l, rule=rule_eq_k4)

        # k5 - Temperature Difference between Tube and Bed
        def rule_eq_k5(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.dThx[i] == b.Ttube[i] - b.Tse[b.feom[i]]
        blk.eq_k5 = Constraint(blk.l, rule=rule_eq_k5)

        # k6 - HX Tube Wall Energy Balance
        def rule_eq_k6(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.ht[i]*b.dThx[i]*b.Cr \
                                    == b.hw*(b.Thx[b.feom[i]]-b.Ttube[i])
        blk.eq_k6 = Constraint(blk.l, rule=rule_eq_k6)
        
        # k7 - HX Fluid Inlet Pressure
        blk.eq_k7 = Constraint(expr = blk.Phx[blk.nfe] == blk.HX_In_P)
        
        # k8 - HX Fluid Inlet Temperature
        blk.eq_k8 = Constraint(expr = blk.Thx[blk.nfe] == blk.HX_In_T)
        
        # k9 - HX Fluid Enthalpy Flow
        def rule_eq_k9(b, i):
            return b.Hhx[i] == b.HX_In_F*b.hx_prop[i].h_mix
        blk.eq_k9 = Constraint(blk.l, rule=rule_eq_k9)

        # ---------------------------------------------------------------------
        # HX Boundary Conditions
        blk.eq_d12 = Constraint(expr = blk.HX_Out_F == blk.HX_In_F)
        blk.eq_d13 = Constraint(expr = blk.HX_Out_T == blk.Thx[0])
        blk.eq_d14 = Constraint(expr = blk.HX_Out_P == blk.Phx[0])
        def rule_d15(b, j):
            return blk.HX_Out_y[j] == blk.HX_In_y[j]
        blk.eq_d15 = Constraint(blk.HXList, rule=rule_d15)
     
    def _make_props(self):
        """
        Create property constraints.
        """
        blk = self
        
        # Create gas property blocks
        def rule_gas_prop_blk_b(b, i):
            block = blk.gas_prop_lib.PropPack(name = 'gprop',
                                              parent = blk,
                                              prop_list = {"V_vap",
                                                           "cp_vap",
                                                           "k_vap",
                                                           "h_vap"})
            return block
        blk.gas_prop_b = Block(blk.l, rule=rule_gas_prop_blk_b)
        
        def rule_gas_prop_blk_c(b, i):
            if i != 0:
                block = blk.gas_prop_lib.PropPack(name = 'gprop',
                                                  parent = blk,
                                                  prop_list = {"h_vap"})
                return block
        blk.gas_prop_c = Block(blk.l, rule=rule_gas_prop_blk_c)
        
        def rule_gas_prop_blk_e(b, i):
            if i != 0:
                block = blk.gas_prop_lib.PropPack(name = 'gprop',
                                                  parent = blk,
                                                  prop_list = {"mu_vap",
                                                               "rho_vap",
                                                               "D_vap",
                                                               "cp_vap",
                                                               "MW_vap",
                                                               "k_vap",
                                                               "h_vap"})
                return block
        blk.gas_prop_e = Block(blk.l, rule=rule_gas_prop_blk_e)
        
        blk.gas_prop_out = blk.gas_prop_lib.PropPack(name = 'gprop',
                                                     parent = blk,
                                                     prop_list = {"h_vap"})

        # Link variables between model and properties
        def rule_gpropb_T(b, i):
            return b.gas_prop_b[i].T == b.Tgb[i]
        blk.eq_gpropb_T = Constraint(blk.l, rule=rule_gpropb_T)
        def rule_gpropb_P(b, i):
            return b.gas_prop_b[i].P == b.P[i]
        blk.eq_gpropb_P = Constraint(blk.l, rule=rule_gpropb_P)
        def rule_gpropb_y(b, i, j):
            return b.gas_prop_b[i].y[j] == b.yb[j,i]
        blk.eq_gpropb_y = Constraint(blk.l, blk.GasList, rule=rule_gpropb_y)
        
        def rule_gpropc_T(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_c[i].T == b.Tgc[i]
        blk.eq_gpropc_T = Constraint(blk.l, rule=rule_gpropc_T)
        def rule_gpropc_P(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_c[i].P == b.P[i]
        blk.eq_gpropc_P = Constraint(blk.l, rule=rule_gpropc_P)
        def rule_gpropc_y(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_c[i].y[j] == b.yc[j,i]
        blk.eq_gpropc_y = Constraint(blk.l, blk.GasList, rule=rule_gpropc_y)

        def rule_gprope_T(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_e[i].T == b.Tge[i]
        blk.eq_gprope_T = Constraint(blk.l, rule=rule_gprope_T)
        def rule_gprope_P(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_e[i].P == b.P[i]
        blk.eq_gprope_P = Constraint(blk.l, rule=rule_gprope_P)
        def rule_gprope_y(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.gas_prop_e[i].y[j] == b.ye[j,i]
        blk.eq_gprope_y = Constraint(blk.l, blk.GasList, rule=rule_gprope_y)
        
        blk.eq_gpropo_T = Constraint(expr = blk.Gas_Out_T \
                                                 == blk.gas_prop_out.T)
        blk.eq_gpropo_P = Constraint(expr = blk.Gas_Out_P \
                                                 == blk.gas_prop_out.P)
        def rule_gpropo_y(b, j):
            return b.Gas_Out_y[j] == b.gas_prop_out.y[j]
        blk.eq_gpropo_y = Constraint(blk.GasList, rule=rule_gpropo_y)

        # ---------------------------------------------------------------------
        # Create solid property blocks
        def rule_sol_prop_blk_c(b, i):
            if i != 0:
                block = blk.sol_prop_lib.PropPack(name = 'sprop',
                                                  parent = blk,
                                                  prop_list = {"v_mf",
                                                               "ap",
                                                               "h_sol",
                                                               "r",
                                                               "rho_sol"})
                return block
        blk.sol_prop_c = Block(blk.l, rule=rule_sol_prop_blk_c)

        def rule_sol_prop_blk_e(b, i):
            if i != 0:
                block = blk.sol_prop_lib.PropPack(name = 'sprop',
                                                  parent = blk)
                return block
        blk.sol_prop_e = Block(blk.l, rule=rule_sol_prop_blk_e)

        def rule_sol_prop_blk_in(b):
            block = blk.sol_prop_lib.PropPack(name = 'sprop',
                                              parent = blk,
                                              prop_list = {"h_sol"})
            return block
        blk.sol_prop_f = Block(rule=rule_sol_prop_blk_in)
        
        # Reaction rate constraints
        def rule_l1(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.rgc[j,i] == b.sol_prop_c[i].rg[j]*b.dl[i]*b.Ax \
                                    * b.delta[i]*b.fcw[i]*(1-b.ed[i])
        blk.eq_l1 = Constraint(blk.l, blk.GasList, rule=rule_l1)

        def rule_l2(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.rsc[j,i] == b.sol_prop_c[i].rs[j]*b.dl[i]*b.Ax \
                                    * b.delta[i]*b.fcw[i]*(1-b.ed[i])
        blk.eq_l2 = Constraint(blk.l, blk.SolidList,rule=rule_l2)

        def rule_l3(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.rge[j,i] == b.sol_prop_e[i].rg[j]*b.dl[i]*b.Ax \
                                    * b.delta_e[i] \
                                    * (1-b.ed[i])
        blk.eq_l3 = Constraint(blk.l, blk.GasList, rule=rule_l3)
        def rule_l4(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.rse[j,i] == b.sol_prop_e[i].rs[j]*b.dl[i]*b.Ax \
                                    * b.delta_e[i] \
                                    * (1-b.ed[i])
        blk.eq_l4 = Constraint(blk.l, blk.SolidList,rule=rule_l4)
        
        # Link variables between model and properties
        def rule_spropc_T(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_c[i].T == b.Tsc[i]
        blk.eq_spropc_T = Constraint(blk.l, rule=rule_spropc_T)
        def rule_spropc_P(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_c[i].P == b.P[i]
        blk.eq_spropc_P = Constraint(blk.l, rule=rule_spropc_P)
        def rule_spropc_y(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_c[i].y[j] == b.yc[j,i]
        blk.eq_spropc_y = Constraint(blk.l, blk.GasList, rule=rule_spropc_y)
        def rule_spropc_x(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_c[i].w[j] == b.xc[j,i]
        blk.eq_spropc_x = Constraint(blk.l, blk.SolidList, rule=rule_spropc_x)

        def rule_sprope_T(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_e[i].T == b.Tse[b.feom[i]]
        blk.eq_sprope_T = Constraint(blk.l, rule=rule_sprope_T)
        def rule_sprope_P(b, i):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_e[i].P == b.P[i]
        blk.eq_sprope_P = Constraint(blk.l, rule=rule_sprope_P)
        def rule_sprope_y(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_e[i].y[j] == b.ye[j,i]
        blk.eq_sprope_y = Constraint(blk.l, blk.GasList, rule=rule_sprope_y)
        def rule_sprope_x(b, i, j):
            if i == 0:
                return Constraint.Skip
            else:
                return b.sol_prop_e[i].w[j] == b.xe[j,b.feom[i]]
        blk.eq_sprope_x = Constraint(blk.l, blk.SolidList, rule=rule_sprope_x)

        def rule_spropf_T(b):
            return b.sol_prop_f.T == b.Solid_In_T
        blk.eq_spropf_T = Constraint(rule=rule_spropf_T)
        def rule_spropf_x(b, j):
            return b.sol_prop_f.w[j] == b.Solid_In_x[j]
        blk.eq_spropf_x = Constraint(blk.SolidList, rule=rule_spropf_x)
        def rule_spropf_P(b):           # NEEDED CURRENTLY FOR SQUARE PROBLEM
            return b.sol_prop_f.P == b.Gas_In_P
        blk.eq_spropf_P = Constraint(rule=rule_spropf_P)
        def rule_spropf_y(b, j):
            return b.sol_prop_f.y[j] == b.Gas_In_y[j]
        blk.eq_spropf_y = Constraint(blk.GasList, rule=rule_spropf_y)

        # ---------------------------------------------------------------------
        # Create HX fluid property blocks
        def rule_hx_prop_blk(b, i):
            block = blk.hx_prop_lib.PropPack(name='hprop',
                                                  parent=blk)
            return block
        blk.hx_prop = Block(blk.l, rule=rule_hx_prop_blk)

        blk.hx_prop_in = blk.hx_prop_lib.PropPack(name='hprop',
                                                  parent=blk)

        # Link variables between model and properties
        def rule_hprop_T(b, i):
            return b.hx_prop[i].T == b.Thx[i]
        blk.eq_hprop_T = Constraint(blk.l, rule=rule_hprop_T)
        def rule_hprop_P(b, i):
            return b.hx_prop[i].P == b.Phx[i]
        blk.eq_hprop_P = Constraint(blk.l, rule=rule_hprop_P)
        def rule_hprop_y(b, i, j):
            return b.hx_prop[i].y[j] == b.HX_In_y[j]
        blk.eq_hprop_y = Constraint(blk.l, blk.HXList, rule=rule_hprop_y)
        
        blk.eq_gprop0_T = Constraint(expr = blk.hx_prop_in.T == blk.HX_In_T)
        blk.eq_gprop0_P = Constraint(expr = blk.hx_prop_in.P == blk.HX_In_P)
        def rule_hprop0_y(b, j):
            return b.hx_prop_in.y[j] == b.HX_In_y[j]
        blk.eq_hprop0_y = Constraint(blk.HXList, rule=rule_hprop0_y)

    def _dae_transform(self):
        # Apply DAE transformations
        blk = self

        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(blk, nfe=blk.nfe, wrt=blk.l, scheme='BACKWARD')
        
        #discretizer = TransformationFactory('dae.collocation')
        #discretizer.apply_to(blk, nfe=blk.nfe, ncp=3, wrt=blk.l,
        #                                          scheme='LAGRANGE-RADAU')
        
        blk.feom = {}
        blk.feop = {}
        for i in range(1,len(blk.l)):
            blk.feom[blk.l[i+1]] = blk.l[i]
            blk.feop[blk.l[i]] = blk.l[i+1]
            
            #if blk.l[i] not in blk.fe_set.keys():
            #    c = blk.l[i]
            #    a = blk.fe_set[ceil(blk.l[i])]
            #    b = blk.fe_set[floor(blk.l[i])]
            #    blk.fe_set[c] =  (a-b)*(c-floor(c)) + b

    def _initialize(blk, Tg=None, Ts=None, Th=None, P=None, x=None, y=None,
                        outlvl=0, optarg=None):
        # Start time for total elapsed time
        total_time_start = time.time()
        
        if optarg == None:
            sopt = {"tol"            : 1e-8,
                    "max_cpu_time"   : 300,
                    "print_level"    : 5}
        else:
            sopt = optarg
        
        if outlvl > 1:
            stee = True
        else:
            stee = False

        if outlvl > 0:
            print("\n")
            print("----------------------------------------------------------")
            print("BFB Reactor Initialisation\n")

        # Set Initial Values of State and Property Variables
        for i in blk.l:
            blk.P[i] = value(blk.Gas_In_P)
            blk.Tgb[i] = value(blk.Gas_In_T)
            blk.Tgc[i] = value(blk.Gas_In_T)
            blk.Tge[i] = value(blk.Gas_In_T)
            blk.Tsc[i] = value(blk.Solid_In_T)
            blk.Tse[i] = value(blk.Solid_In_T)
            for j in blk.GasList:
                blk.yb[j,i] = value(blk.Gas_In_y[j])
                blk.yc[j,i] = value(blk.Gas_In_y[j])
                blk.ye[j,i] = value(blk.Gas_In_y[j])
            for j in blk.SolidList:
                blk.xc[j,i] = value(blk.Solid_In_x[j])
                blk.xe[j,i] = value(blk.Solid_In_x[j])

            blk.Thx[i] = value(blk.HX_In_T)
            blk.Phx[i] = value(blk.HX_In_P)


        # Property Initial Values
        for i in blk.l:
            blk.gas_prop_b[i].V_vap = value(blk.R)*value(blk.Tgb[i]) \
                                    / value(blk.P[i])
            for j in blk.GasList:
                blk.cb[j,i] = value(blk.Gas_In_y[j]) \
                                / value(blk.gas_prop_b[i].V_vap)
                blk.cc[j,i] = value(blk.cb[j,i])
                blk.ce[j,i] = value(blk.cb[j,i])

        # ---------------------------------------------------------------------
        # Deactivate constraints
        # All property blocks active

        # All group A constraints active

        blk.eq_b1.deactivate()
        blk.eq_b2.deactivate()
        blk.eq_b3.deactivate()
        blk.eq_b4.deactivate()
        blk.eq_b5.deactivate()
        blk.eq_b6.deactivate()
        blk.eq_b7.deactivate()
        blk.eq_b8.deactivate()
        blk.eq_b9.deactivate()
        blk.eq_b10.deactivate()
        blk.eq_b11.deactivate()
        blk.eq_b12.deactivate()
        blk.eq_b13.deactivate()
        blk.eq_b14.deactivate()
        blk.eq_b15.deactivate()
        blk.eq_b16.deactivate()
        blk.eq_b17.deactivate()
        blk.eq_b18.deactivate()
        blk.eq_b19.deactivate()
        blk.eq_b20.deactivate()

        blk.eq_c1.deactivate()
        # c2 and c3 active
        blk.eq_c4.deactivate()
        blk.eq_c5.deactivate()
        blk.eq_c6.deactivate()
        blk.eq_c7.deactivate()
        blk.eq_c8.deactivate()
        blk.eq_c9.deactivate()
        blk.eq_c10.deactivate()
        blk.eq_c11.deactivate()
        blk.eq_c12.deactivate()
        blk.eq_c13.deactivate()
        blk.eq_c14.deactivate()
        blk.eq_c15.deactivate()

        # d1 - d4 active
        blk.eq_d5.deactivate()
        blk.eq_d6.deactivate()
        blk.eq_d7.deactivate()
        blk.eq_d8.deactivate()
        # d9 - d11 active
        blk.eq_d12.deactivate()
        blk.eq_d13.deactivate()
        blk.eq_d14.deactivate()
        blk.eq_d15.deactivate()

        blk.eq_e1.deactivate()
        blk.eq_e2.deactivate()
        blk.eq_e3.deactivate()
        blk.eq_e4.deactivate()
        blk.eq_e5.deactivate()
        blk.eq_e6.deactivate()

        blk.eq_f1.deactivate()
        blk.eq_f2.deactivate()
        blk.eq_f3.deactivate()
        blk.eq_f4.deactivate()
        blk.eq_f5.deactivate()
        blk.eq_f6.deactivate()
        blk.eq_f7.deactivate()
        blk.eq_f8.deactivate()
        blk.eq_f9.deactivate()
        blk.eq_f10.deactivate()

        blk.eq_g1.deactivate()

        blk.eq_h1.deactivate()
        blk.eq_h2.deactivate()
        blk.eq_h3.deactivate()
        blk.eq_h4.deactivate()
        blk.eq_h5.deactivate()
        blk.eq_h6.deactivate()
        blk.eq_h7.deactivate()
        blk.eq_h8.deactivate()
        blk.eq_h9.deactivate()
        blk.eq_h10.deactivate()
        blk.eq_h11.deactivate()
        blk.eq_h12.deactivate()

        blk.eq_i1.deactivate()
        blk.eq_i2.deactivate()
        blk.eq_i3.deactivate()
        blk.eq_i4.deactivate()
        blk.eq_i5.deactivate()
        blk.eq_i6.deactivate()

        blk.eq_j1.deactivate()
        blk.eq_j2.deactivate()
        blk.eq_j3.deactivate()
        blk.eq_j4.deactivate()
        blk.eq_j5.deactivate()
        blk.eq_j6.deactivate()
        blk.eq_j7.deactivate()
        blk.eq_j8.deactivate()

        blk.eq_k1.deactivate()
        blk.eq_k2.deactivate()
        blk.eq_k3.deactivate()
        blk.eq_k4.deactivate()
        blk.eq_k5.deactivate()
        blk.eq_k6.deactivate()
        blk.eq_k7.deactivate()
        blk.eq_k8.deactivate()
        # k9 active

        blk.eq_l1.deactivate()
        blk.eq_l2.deactivate()
        blk.eq_l3.deactivate()
        blk.eq_l4.deactivate()

        # ---------------------------------------------------------------------
        # 1st Initialisation Step - Properties Initialisation
        # Initialise Variables
        for i in blk.l:
            if i != 0:
                for j in blk.GasList:
                    blk.rgc[j,i] = 0.0
                    blk.rge[j,i] = 0.0
                for j in blk.SolidList:
                    blk.rsc[j,i] = 0.0
                    blk.rse[j,i] = 0.0

        # Fix Variables
        for i in blk.l:
            blk.P[i].fix()
            if i != 0:
                blk.Tgb[i].fix()
                blk.Tgc[i].fix()
                blk.Tge[i].fix()
                for j in blk.GasList:
                    blk.yb[j,i].fix()
                    blk.yc[j,i].fix()
                    blk.ye[j,i].fix()
                    blk.cc[j,i].fix()
                    blk.ce[j,i].fix()
                blk.sol_prop_c[i].r.fix(0.0)
                blk.sol_prop_e[i].r.fix(0.0)
                blk.Qhx[i].fix(0.0)
        blk.Tsc.fix()
        blk.Tse.fix()
        blk.Thx.fix()
        blk.Phx.fix()
        
        blk.Gas_Out_F.fix(value(blk.Gas_In_F))
        blk.Gas_Out_T.fix(value(blk.Gas_In_T))
        blk.Gas_Out_P.fix(value(blk.Gas_In_P))
        for j in blk.GasList:
            blk.Gas_Out_y[j].fix(value(blk.Gas_In_y[j]))

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  1, Initialise Properties,                     Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 2nd Initialisation Step - Hydrodynamics
        # Initialise variables
        for i in blk.l:
            blk.vg[i] = value(blk.Gas_In_F)*(value(blk.R)*value(blk.Gas_In_T) \
                                / value(blk.Gas_In_P)) \
                            / (0.25*value(blk.pi)*value(blk.Dt)**2)
            if i != 0:
                blk.Gb[i] = value(blk.Gb[0])
                blk.Ge[i] = value(blk.Ge[0])
                blk.ve[i] = (1/3)*(value(blk.vg[i]) \
                                    -value(blk.sol_prop_e[i].v_mf)) \
                                + value(blk.sol_prop_e[i].v_mf)
                blk.dbm[i] = 2.59*(value(blk.gc)**(-0.2)) \
                                * ((value(blk.vg[i])-value(blk.ve[i])) \
                                * value(blk.Ax))**0.4
                blk.g1[i] = 2.56E-2*sqrt(value(blk.Dt)/value(blk.gc)) \
                                / value(blk.sol_prop_e[i].v_mf)
        for i in blk.l:
            if i == 0:
                blk.db[0] = 1.38*(value(blk.gc)**(-0.2)) \
                            * ((value(blk.vg[blk.feop[i]]) \
                        - value(blk.ve[blk.feop[i]]))*value(blk.Ao))**0.4
            else:
                db_a = 1+0.3*value(blk.dl[i])/value(blk.Dt)
                db_b = 0.3*value(blk.dl[i])*value(blk.g1[i]) \
                            / sqrt(value(blk.Dt))
                db_c = -0.3*value(blk.dl[i])*value(blk.dbm[i])/value(blk.Dt) \
                        - value(blk.db[blk.feom[i]])
            
                blk.db[i] = ((-db_b + sqrt(db_b**2 - 4*db_a*db_c))/(2*db_a))**2
            if i != 0:
                blk.vbr[i] = 0.711*sqrt(value(blk.gc)*value(blk.db[i]))
                if blk.vb_method == 'Werther A':
                    blk.vb[i] = 1.55*((value(blk.vg[i]) \
                            - value(blk.sol_prop_e[i].v_mf)) \
                            + 14.1*(value(blk.db[i])+0.005)) \
                                * (value(blk.Dte)**0.32) + value(blk.vbr[i])
                elif blk.vb_method == 'Werther B':
                    blk.vb[i] = 1.6*((value(blk.vg[i]) \
                            - value(blk.sol_prop_e[i].v_mf)) \
                            + 1.13*(value(blk.db[i])**0.5)) \
                                * (value(blk.Dte)**1.35) + value(blk.vbr[i])
                else:
                    blk.vb[i] = value(blk.vg[i]) + value(blk.vbr[i]) \
                              - value(blk.sol_prop_e[i].v_mf)
                blk.delta[i] = value(blk.Gas_In_F) \
                            * value(blk.gas_prop_b[i].V_vap) \
                            / (value(blk.vb[i])*value(blk.Ax))
                blk.fc[i] = -6*(value(blk.sol_prop_e[i].v_mf) \
                            / value(blk.sol_prop_e[i].e_mf)) \
                            / ((1-3/value(blk.fc_max)) \
                            * (value(blk.sol_prop_e[i].v_mf) \
                            / value(blk.sol_prop_e[i].e_mf)) \
                            - value(blk.vbr[i]) \
                            - (((1+3/value(blk.fc_max)) \
                                * (value(blk.sol_prop_e[i].v_mf) \
                                    / value(blk.sol_prop_e[i].e_mf)) \
                                - value(blk.vbr[i]))**2 \
                                + value(blk.eps))**0.5)
                blk.fcw[i] = value(blk.fc[i]) + value(blk.fw)
                blk.ed[i] = 1-0.958773*2.05*(1-value(blk.sol_prop_e[i].e_mf)) \
                            * (value(blk.sol_prop_e[i].dp)**0.1) \
                            * (value(blk.gc)**0.118) \
                            / (2.54*value(blk.gas_prop_e[i].mu_vap)**0.066) \
                            * value(blk.Lb)**0.043
                blk.e[i] = 1- (1 - value(blk.ed[i]))*(1-value(blk.delta[i]))
                blk.delta_e[i] = 1 - (value(blk.fcw[i])+1)*value(blk.delta[i])

        for i in blk.l:
            if i == 0:
                blk.Jc[0] = value(blk.fw)*value(blk.delta[blk.feop[0]]) \
                            * value(blk.sol_prop_e[blk.feop[0]].rho_sol) \
                            * (1-value(blk.ed[blk.feop[0]])) \
                            * value(blk.vb[blk.feop[0]])
            else:
                blk.Jc[i] = value(blk.fw)*value(blk.delta[i]) \
                            * value(blk.sol_prop_e[i].rho_sol) \
                            * (1-value(blk.ed[i]))*value(blk.vb[i])
            blk.Je[i] = value(blk.Jc[i])
            for j in blk.SolidList:
                blk.Jcc[j,i] = value(blk.Jc[i])*value(blk.xc[j,i])
                blk.Jec[j,i] = value(blk.Je[i])*value(blk.xe[j,i])
            blk.Jch[i] = value(blk.Jc[i])*value(blk.sol_prop_f.h_sol)
            blk.Jeh[i] = value(blk.Je[i])*value(blk.sol_prop_f.h_sol)
            if i != 0:
                blk.Ksbulk[i] = 0.0
                for j in blk.SolidList:
                    blk.Ksbulk_c[j,i] = 0.0
            if i != 0:
                blk.cct[i] = 1/value(blk.gas_prop_b[i].V_vap)
                blk.cet[i] = 1/value(blk.gas_prop_b[i].V_vap)
                for j in blk.GasList:
                    blk.cb[j,i] = value(blk.Gas_In_y[j]) \
                                / value(blk.gas_prop_b[i].V_vap)
                    blk.cc[j,i] = value(blk.cb[j,i])
                    blk.ce[j,i] = value(blk.cb[j,i])
            if i != 0:
                blk.Kcebs[i] = 3*(1-value(blk.ed[i]))*value(blk.ve[i]) \
                                / ((1-value(blk.delta[i]))*value(blk.ed[i]) \
                                * value(blk.db[i]))
                blk.Hbc[i] = 1.32*4.5*value(blk.sol_prop_c[i].v_mf) \
                            * value(blk.gas_prop_b[i].cp_vap) \
                            / (value(blk.db[i]) \
                            * value(blk.gas_prop_b[i].V_vap)) \
                        + 5.85*sqrt(value(blk.gas_prop_b[i].k_vap) \
                            * value(blk.gas_prop_b[i].cp_vap) \
                            / value(blk.gas_prop_b[i].V_vap)) \
                            * (value(blk.gc)**0.25)/value(blk.db[i])**(5/4)
                blk.Hce[i] = 6.78*sqrt(value(blk.ed[i])*value(blk.vbr[i]) \
                        * value(blk.gas_prop_e[i].k_vap)*value(blk.cct[i]) \
                        * value(blk.gas_prop_e[i].cp_vap)) \
                        / value(blk.db[i])**(3/2)
                blk.Hgbulk[i] == (6*value(blk.Kd)*value(blk.delta[i]) \
                            * value(blk.dl[i])*value(blk.Ax) \
                            * (value(blk.cet[i]) \
                        - (1/value(blk.gas_prop_b[i].V_vap))) \
                            / value(blk.db[i])) \
                            * value(blk.gas_prop_b[i].cp_vap) \
                            * value(blk.Tgb[i])
                blk.Hsbulk[i] = value(blk.Ksbulk[i]) \
                            * value(blk.sol_prop_e[i].h_sol)
                for j in blk.GasList:
                    blk.Kbc[j,i] = 1.32*4.5*(value(blk.sol_prop_c[i].v_mf) \
                                / value(blk.db[i])) \
                                + 5.85*((value(blk.gc)**0.25) \
                                * (value(blk.gas_prop_e[i].D_vap[j]))**0.5 \
                                / (value(blk.db[i])**(5/4)))
                    blk.Kce[j,i] = 6.77*sqrt(value(blk.ed[i]) \
                                * value(blk.gas_prop_e[i].D_vap[j]) \
                                * value(blk.vbr[i])/value(blk.db[i])**3)
                for j in blk.SolidList:
                    blk.Kcebs_c[j,i] = value(blk.dl[i])*value(blk.Ax) \
                                * value(blk.delta[i]) \
                                * value(blk.sol_prop_e[i].rho_sol) \
                                * value(blk.Kcebs[i])*(value(blk.xc[j,i]) \
                                - value(blk.xe[j,i]))

        # Fix variables
        blk.vg.fix()
        blk.Jcc.fix()
        blk.Jec.fix()
        blk.Jch.fix()
        blk.Jeh.fix()
        for i in blk.l:
            if i > 0:
                blk.delta[i].fix()
                blk.Gb[i].fix()
                blk.Ge[i].fix()
                blk.cct[i].fix()
                blk.cet[i].fix()

        # Activate constraints
        blk.eq_b7.activate()
        blk.eq_b8.activate()

        blk.eq_c6.activate()
        blk.eq_c7.activate()
        blk.eq_c8.activate()
        blk.eq_c9.activate()
        blk.eq_c10.activate()

        blk.eq_f1.activate()
        blk.eq_f2.activate()
        blk.eq_f3.activate()
        blk.eq_f4.activate()
        blk.eq_f5.activate()
        blk.eq_f6.activate()
        blk.eq_f7.activate()
        blk.eq_f8.activate()
        blk.eq_f9.activate()
        blk.eq_f10.activate()

        blk.eq_g1.activate()

        blk.eq_h1.activate()
        blk.eq_h2.activate()
        blk.eq_h3.activate()
        blk.eq_h4.activate()
        blk.eq_h5.activate()
        blk.eq_h6.activate()
        blk.eq_h7.activate()
        blk.eq_h8.activate()
        blk.eq_h9.activate()
        blk.eq_h10.activate()
        blk.eq_h11.activate()
        blk.eq_h12.activate()

        blk.eq_i1.activate()
        blk.eq_i2.activate()
        blk.eq_i3.activate()
        blk.eq_i4.activate()
        blk.eq_i5.activate()
        blk.eq_i6.activate()

        blk.eq_l1.activate()
        blk.eq_l2.activate()
        blk.eq_l3.activate()
        blk.eq_l4.activate()

        # Unfix variables
        blk.P.unfix()

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  2, Hydrodynamics,                             Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 3rd Initialisation Step - Mass Balances I

        # Activate constraints
        blk.eq_b1.activate()
        blk.eq_b2.activate()
        blk.eq_b3.activate()
        blk.eq_b4.activate()
        blk.eq_b5.activate()
        blk.eq_b6.activate()
        blk.eq_b9.activate()
        blk.eq_b10.activate()
        #blk.eq_b13.activate()
        #blk.eq_b14.activate()
        blk.eq_b17.activate()
        blk.eq_b18.activate()

        blk.eq_c1.activate()
        blk.eq_c4.activate()
        blk.eq_c5.activate()
        blk.eq_c11.activate()
        blk.eq_c12.activate()
        blk.eq_c13.activate()

        blk.eq_e1.activate()
        blk.eq_e2.activate()
        blk.eq_e3.activate()
        blk.eq_e4.activate()
        blk.eq_e5.activate()
        blk.eq_e6.activate()

        # Unfix variables
        blk.delta.unfix()
        blk.Gb.unfix()
        blk.Ge.unfix()
        blk.vg.unfix()
        blk.cct.unfix()
        blk.cet.unfix()
        blk.yb.unfix()
        blk.yc.unfix()
        blk.ye.unfix()
        blk.cc.unfix()
        blk.ce.unfix()
        blk.Tgb.unfix()
        blk.Tgc.unfix()
        blk.Tge.unfix()
        blk.Jcc.unfix()
        blk.Jec.unfix()
        for i in blk.l:
            if i > 0:
                blk.sol_prop_c[i].r.unfix()
                blk.sol_prop_e[i].r.unfix()
                
        if blk.s_inlet == 'Top':
            blk.eq_b13.activate()
            for j in blk.SolidList:
                blk.Jec[j,blk.feom[blk.nfe]].fix(
                            value(blk.Je[blk.feom[blk.nfe]]) \
                            * value(blk.xe[j,blk.feom[blk.nfe]]))
        else:
            blk.eq_b14.activate()
            for j in blk.SolidList:
                blk.Jec[j,0].fix(value(blk.Je[0])*value(blk.xe[j,0]))

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  3, Mass Balances I,                           Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 4th Initialisation Step - Mass Balances II
        # Unfix variables
        blk.Jec.unfix()

        # Activate constraints
        blk.eq_b13.activate()
        blk.eq_b14.activate()

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  4, Mass Balances II,                          Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 5th Initialisation Step - Energy Balances I

        # Activate constraints
        blk.eq_b11.activate()
        blk.eq_b12.activate()
        blk.eq_b19.activate()
        blk.eq_b20.activate()
        
        blk.eq_c14.activate()
        blk.eq_c15.activate()

        # Unfix variables
        blk.Jch.unfix()
        blk.Jeh.unfix()
        blk.Tsc.unfix()
        blk.Tse.unfix()

        if blk.s_inlet == 'Top':
            blk.eq_b15.activate()
            blk.Jeh[blk.feom[blk.nfe]].fix(value(blk.Je[blk.feom[blk.nfe]]) \
                           * value(blk.sol_prop_c[blk.nfe].h_sol))
        else:
            blk.eq_b16.activate()
            blk.Jeh[0].fix(value(blk.Je[0])*value(blk.sol_prop_e[1].h_sol))

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  5, Energy Balances I,                         Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 6th Initialisation Step - Energy Balances II
        # Unfix variables
        blk.Jeh.unfix()

        # Activate constraints
        blk.eq_b15.activate()
        blk.eq_b16.activate()

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  6, Energy Balances II,                        Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 7th Initialisation Step - HX Tubes
        # Initialise variables
        for i in blk.l:
            if i > 0:
                blk.kpa[i] = (3.58-2.5*value(blk.ed[i])) \
                                * value(blk.gas_prop_e[i].k_vap) \
                                * ((value(blk.sol_prop_e[i].k_sol) \
                                    / value(blk.gas_prop_e[i].k_vap)) \
                                ** (0.46-0.46*value(blk.ed[i])))
                blk.fn[i] = value(blk.vg[i])/value(blk.sol_prop_e[i].v_mf)
                blk.tau[i] = 0.44*((value(blk.sol_prop_e[i].dp)*value(blk.gc) \
                                / ((value(blk.sol_prop_e[i].v_mf)**2) \
                                * ((value(blk.fn[i]) \
                                    - value(blk.ah))**2)))**0.14) \
                                * ((value(blk.sol_prop_e[i].dp) \
                                    / value(blk.dx))**0.225)
                blk.fb[i] = 0.33*(((value(blk.sol_prop_e[i].v_mf)**2) \
                                * ((value(blk.fn[i])-value(blk.ah))**2) \
                                / (value(blk.sol_prop_e[i].dp) \
                                * value(blk.gc)))**0.14)
                blk.hd[i] = 1000
                blk.hl[i] = 10
                blk.ht[i] = 100

        # Activate constraints
        blk.eq_d12.activate()
        blk.eq_d13.activate()
        blk.eq_d14.activate()
        blk.eq_d15.activate()

        blk.eq_j1.activate()
        blk.eq_j2.activate()
        blk.eq_j3.activate()
        blk.eq_j4.activate()
        blk.eq_j5.activate()
        blk.eq_j6.activate()
        blk.eq_j7.activate()
        blk.eq_j8.activate()

        blk.eq_k1.activate()
        blk.eq_k2.activate()
        blk.eq_k3.activate()
        blk.eq_k4.activate()
        blk.eq_k5.activate()
        blk.eq_k6.activate()
        blk.eq_k7.activate()
        blk.eq_k8.activate()

        # Unfix variables
        blk.Qhx.unfix()
        blk.Thx.unfix()
        blk.Phx.unfix()

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  7, HX Tubes,                                  Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        # 8th Initialisation Step - Outlets
        # Activate constraints
        blk.eq_d5.activate()
        blk.eq_d6.activate()
        blk.eq_d7.activate()
        blk.eq_d8.activate()

        # Unfix variables
        blk.Gas_Out_F.unfix()
        blk.Gas_Out_T.unfix()
        blk.Gas_Out_P.unfix()
        blk.Gas_Out_y.unfix()

        ts = time.time()
        # Solve and print progress
        if outlvl > 0:
            print("Step  8, Outlets,                                   Time: ",
                  end="")
        results = blk.solve(options=sopt,tee=stee)
        if outlvl > 0:
            print("{}, {}".format(model_util.hhmmss(time.time() - ts),
                    results.solver.message))

        # ---------------------------------------------------------------------
        if outlvl > 0:        
            print("\n\nTotal Initialisation Time: {0}".format(
                model_util.hhmmss(time.time()-total_time_start)))
            print("BFB Initialization Complete")
            print("----------------------------------------------------------")

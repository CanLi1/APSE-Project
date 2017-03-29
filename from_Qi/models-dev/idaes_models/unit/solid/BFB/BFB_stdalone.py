# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 13:26:46 2016

@author: alee

PETSc version, SI units
"""

# Import Necessary Packages
from __future__ import division
import time

# Model start time for duration calculations
stime = time.time()
print('Importing packages')

from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus
from pyomo.dae import *

print('Import Complete                 (','%.5f'%(time.time()-stime),'s)')
splitt = time.time()
print('Beginning Model Construction')

# =============================================================================
# Declare Model
model = ConcreteModel()

metadata = {'Model':'BFB',
            'Creator':'A Lee',
            'Date':'2016-09-13',
            'Version':'7.0.0.1',
            'Pyomo':'4.4.1',
            'IDAES Standards':'0.3'}

# -----------------------------------------------------------------------------
# Declare Component Lists
GasList = model.GasList = ['CO2','H2O','N2']
SolidList = model.SolidList = ['H2O','Bic','Car']
HXList = model.HXList = ['H2O']

# -----------------------------------------------------------------------------
# Declare Parameters

pi = model.pi = Param(default=2*acos(0), doc='pi')
R = model.R   = Param(default=8.314472,
                      doc='Universal Gas Constant [J/molo.K]')
gc = model.gc = Param(default=9.81,
                      doc='Gravitational Acceleration Constant [m^2/s]')
SmoothIF_eps = modelSmoothIF_eps = Param(default=1e-4,
                      doc='Smoothing Factor for Smooth IF Statements')
Tref = model.Tref = Param(default=298.15,
                      doc='THermodynamic Reference Temperature [K]')

ah = model.ah = Param(within=PositiveReals, mutable=True, default=0.8,
                      doc='Emprical Factor for HX Tube Model')
Cr = model.Cr = Param(within=PositiveReals, mutable=True, default=1,
                      doc='Average Correction Factor for HX Tube Model')
fw = model.fw = Param(within=PositiveReals, mutable=True, default=0.2,
                      doc='Bubble to Wake Volume Ratio')
hw = model.hw = Param(within=PositiveReals, mutable=True, default=1500,
                    doc='Heat Transfer Coefficient of HX Tube Walls [W/m^2.K]')
Kd = model.Kd = Param(within=PositiveReals, mutable=True, default=1,
                      doc='Bulk Gas Permeation Coefficient [m/s]')

# -----------------------------------------------------------------------------
# Declare Distribution Domain
ND = model.ND = Param(default=60, doc='Number of Discrete Element in Model')

l = model.l = RangeSet(1, ND, doc='Set of Discrete Elements')

dl = model.dl = Var(domain=PositiveReals, doc='Size of Discrete Elements [m]')

# -----------------------------------------------------------------------------
# Declare Variables

# Vessel Dimensions
Ax = model.Ax = Var(domain=PositiveReals,
                    doc='Cross-Sectional Area of Bed [m^2]')
Areact = model.Areact = Var(domain=PositiveReals,
                    doc='Cross-Sectional Area of Reactor Vessel [m^2]')
Dt = model.Dt = Var(domain=PositiveReals,
                    doc='Diameter of Reactor Vessel [m]')
Dte = model.Dte = Var(domain=PositiveReals,
                    doc='Hydraulic Diameter of Bed [m]')
Lb = model.Lb = Var(domain=PositiveReals,
                    doc='Depth of Bed [m]')

# Distributor Design
Ao = model.Ao = Var(domain=PositiveReals,
                    doc='Distributor Plate Area per Orifice [m^2/orifice]')
nor = model.nor = Var(domain=PositiveReals,
                      doc='Distributor Plate Orifices per Area [orifices/m^2]')

# Heat Exchanger Dimensions
Ahx = model.Ahx = Var(domain=PositiveReals,
                      doc='Total Area of Heat Exchanger Surfaces [m^2]')
dx = model.dx = Var(domain=PositiveReals,
                    doc='Diamete of Heat Exchanger Tubes [m]')                    
lhx = model.lhx = Var(domain=PositiveReals,
                      doc='Heat Exchanger Tube Spacing (Pitch-Diameter) [m]')
lp = model.lp = Var(domain=PositiveReals,doc='Heat Exchanger Tube Pitch [m]')
Nx = model.Nx = Var(domain=PositiveReals,
                    doc='Number of Heat Exchanger Tubes [-]')

# Gas Inlet Conditions
Gas_In_F = model.Gas_In_F = Var(domain=PositiveReals,
                    doc='Molar Flowrate of Gas at Gas Inlet [mol/s]')
Gas_In_P = model.Gas_In_P = Var(domain=PositiveReals,
                    doc='Pressure of Gas at Gas Inlet [Pa]')
Gas_In_T = model.Gas_In_T = Var(domain=PositiveReals,
                    doc='Temperature of Gas at Gas Inlet [K]')
Gas_In_y = model.Gas_In_y = Var(GasList, domain=PositiveReals,
                    doc='Mole Fractions of Gas Species at Gas Inlet [mol/mol]')

# Solids Inlet Conditions
Solid_In_F = model.Solid_In_F = Var(domain=PositiveReals,
                    doc='Mass Flowrate of Solids at Solid Inlet [kg/s]')
Solid_In_T = model.Solid_In_T = Var(domain=PositiveReals,
                    doc='Temperature of Solids at Solids Inlet [k]')
Solid_In_w = model.Solid_In_w = Var(SolidList, domain=PositiveReals,
                    doc='Loading of Solids at Solids Inlet [mol/kg]')

# HX Fluid Inlet Conditions
HX_In_F = model.HX_In_F = Var(domain=PositiveReals,
                    doc='Heat Exchanger Fluid Flowrate at Inlet [mol/s]')
HX_In_T = model.HX_In_T = Var(domain=PositiveReals,
                    doc='Heat Exchanger Fluid Temperature at Inlet [K]')
HX_In_P = model.HX_In_P = Var(domain=PositiveReals,
                    doc='Heat Exchanger Fluid Pressure at Inlet [Pa]')
HX_In_y = model.HX_In_y = Var(HXList, domain=PositiveReals,
                    doc='Heat Exchanger Fluid Mole Fractions [mol/mol]')

# Material Flows
Gb = model.Gb = Var(l, domain=PositiveReals,
                    doc='Bubble Region Molar Gas FLowrate [mol/s]')
Jc = model.Jc = Var(l, domain=PositiveReals,
                    doc='Cloud-Wake Region Solids Flux [kg/m^2.s]')
Je = model.Je = Var(l, domain=Reals,
                    doc='Emulsion Region Solids Flux (downwards) [kg/m^2.s]')

# Enthalpy Flows
Ghb = model.Ghb = Var(l,doc='Bubble Region Gas Enthalpy Flowrate [J/s]')
Jhc = model.Jhc = Var(l,doc='Cloud-Wake Region Solids Enthalpy Flux [J/m^2.s]')
Jhe = model.Jhe = Var(l,doc='Emulsion Region Solids Enthalpy Flux [J/m^2.s]')

# Temperatures and Pressures
P = model.P = Var(l, domain=PositiveReals, doc='Pressure [Pa]')
Tgb = model.Tgb = Var(l, domain=PositiveReals,
                      doc='Bubble Region Gas Temperature [K]')
Tgc = model.Tgc = Var(l, domain=PositiveReals,
                      doc='Cloud-Wake Region Gas Temperature [K]')
Tge = model.Tge = Var(l, domain=PositiveReals,
                      doc='Emulsion Region Gas Temperature [K]')
Tsc = model.Tsc = Var(l, domain=PositiveReals,
                      doc='Cloud-Wake Region Solids Temperature [K]')
Tse = model.Tse = Var(l, domain=PositiveReals,
                      doc='Emulsion Region Solids Temperature [K]')

# Gas and Solid Compositions
yb = model.yb = Var(GasList, l, domain=PositiveReals,
                    doc='Bubble Region Gas Mole Fractions [mol/mol]')
yc = model.yc = Var(GasList, l, domain=PositiveReals,
                    doc='Cloud-Wake Region Gas Mole Fractions [mol/mol]')
ye = model.ye = Var(GasList, l, domain=PositiveReals,
                    doc='Emulsion Region Gas Mole Fractions [mol/mol]')
wc = model.wc = Var(SolidList, l, domain=PositiveReals,
                    doc='Cloud-Wake Region Solids Loading [mol/kg]')
we = model.we = Var(SolidList, l, domain=PositiveReals,
                    doc='Emulsion Region Solids Loading [mol/kg]')
                    
# Gas Phase Concentrations
cb = model.cb = Var(GasList, l, domain=PositiveReals,
                doc='Bubble Region Gas Component Concentrations [mol/m^3]')
cc = model.cc = Var(GasList, l, domain=PositiveReals,
                doc='Cloud-Wake Region Gas Component Concentrations [mol/m^3]')
ce = model.ce = Var(GasList, l, domain=PositiveReals,
                doc='Emulsion Region Gas Component Concentrations [mol/m^3]')
cct = model.cct = Var(l, domain=PositiveReals,
                doc='Cloud-Wake Region Total Gas Concentration [mol/m^3]')
cet = model.cet = Var(l, domain=PositiveReals,
                doc='Emulsion Region Total Gas Concentration [mol/m^3]')

# Velocities
vb = model.vb = Var(l, domain=PositiveReals,
                    doc='Bubble Velocity [m/s]')
vg = model.vg = Var(l, domain=PositiveReals,
                    doc='Superficial Gas Veleocity [m/s]')
us = model.us = Var(l, doc='Emulsion Region Solids Velocity (downwards) [m/s]')

# Fluidisation Properties
Ar = model.Ar = Var(l, domain=PositiveReals,
                    doc='Archimedes Number [-]')
vmf = model.vmf = Var(domain=PositiveReals,
                    doc='Solids Minimum Fluidisation Velocity [m/s]')

# Bubble Dimensions and Hydrodynamics
db = model.db = Var(l, domain=PositiveReals,
                    doc='Average Bubble Diameter [m]')
db0 = model.db0 = Var(domain=PositiveReals,
                    doc='Average Bubble Diameter above Distributor [m]')
dbm = model.dbm = Var(l, domain=PositiveReals,
                    doc='Maximum Theoretical Bubble Diameter [m]')
delta = model.delta = Var(l, domain=PositiveReals,
                    doc='bubble Region Volume Fraction [m^3/m^3]')
e = model.e = Var(l, domain=PositiveReals,
                    doc='Cross-sEctional Average Voidage Fraction [m^3/m^3]')
ed = model.ed = Var(l, domain=PositiveReals,
                    doc='Emulsion Region Voidage Fraction [m^3/m^3]')
fc = model.fc = Var(l, domain=PositiveReals,
                    doc='Cloud to Bubble Region Volume Ratio [m^3/m^3]')
fcw = model.fcw = Var(l, domain=PositiveReals,
                    doc='Cloud-Wake to Bubble Region Volume Ratio [m^3/m^3]')
g1 = model.g1 = Var(domain=PositiveReals,
                    doc='Bubble Growth Coefficient [-]')
vbr = model.vbr = Var(l, domain=PositiveReals,
                    doc='Bubble Rise Velocity [m/s]')
ve = model.ve = Var(l, doc='Emulsion Region Gas Veloicty [m/s]')

# Bulk Material Transfer
Kgbulk_c = model.Kgbulk_c = Var(GasList, l,
                    doc='Gas Phase Component Bulk Transfer Rate [mol/s]')
Ksbulk = model.Ksbulk = Var(l, doc='Solid Phase Bulk Transfer Rate [kg/s]')
Ksbulk_c = model.Ksbulk_c = Var(SolidList, l,
                    doc='Adsorbed Species Bulk Transfer Rate [mol/s]')
Hgbulk = model.Hgbulk = Var(l,
                    doc='Gas Phase Bulk Enthalpy Trasnfer Rate [J/s]')
Hsbulk = model.Hsbulk = Var(l,
                    doc='Solid Phase Bulk Enthalpy Trasnfer Rate [J/s]')

# Heat and Mass Transfer Coefficients
Kbc = model.Kbc = Var(GasList, l, domain=PositiveReals,
    doc='Bubble to Cloud-Wake Gas Mass Transfer Coeffiicent [1/s]')
Kce = model.Kce = Var(GasList, l, domain=PositiveReals,
    doc='Cloud-Wake to Emulsion Gas Mass Transfer Coefficient [1/s]')
Kcebs = model.Kcebs = Var(l,
    doc='Cloud-Wake to Emulsion Solid Mass Transfer Coefficient [1/s]')
Kcebs_c = model.Kcebs_c = Var(SolidList,l, domain=Reals,
    doc='Cloud-Wake to Emulsion Adsorbed Species Mass Transfer Rate [mol/s]')
Hbc = model.Hbc = Var(l, domain=PositiveReals,
    doc='Bubble to Cloud-Wake Gas Energy Transfer Coefficient [J/m^3.K.s]')
Hce = model.Hce = Var(l, domain=PositiveReals,
    doc='Cloud-Wake to Emulsion Gas Energy Transfer Coefficient [J/m^3.K.s]')
hp = model.hp = Var(l, domain=PositiveReals,
    doc='Gas to Solid Energy Convective Heat Transfer Coefficient [J/m^2.K.s]')

# HX Fluid State
Phx = model.Phx = Var(l, domain=PositiveReals,
                    doc='Heat Exchanger Fluid Pressure [Pa]')
Thx = model.Thx = Var(l, domain=PositiveReals,
                    doc='Heat Exchanger Fluid Temperature [K]')
Ttube = model.Ttube = Var(l, domain=PositiveReals,
                    doc='Heat Exchanger Tube Temperature [K]')
dThx = model.dThx = Var(l,
                    doc='Temperature Difference between HX Tubes and Bed [K]')

# HX Tube Heat Transfer Coefficients
fb = model.fb = Var(l, domain=PositiveReals,
    doc='Fraction of Time HX Tubes Contact DEnse Packets [-]')
fn = model.fn = Var(l, domain=PositiveReals,
    doc='Fluidisation Number [-]') 
hd = model.hd = Var(l, domain=PositiveReals,
    doc='Convective Heat Transfer Coeffiicent of Dense Packets [J/m^2.K.s]')
hl = model.hl = Var(l, domain=PositiveReals,
    doc='Convective Heat Transfer Coeffiicent of Gas Bubbles [J/m^2.K.s]') 
ht = model.ht = Var(l, domain=PositiveReals,
    doc='Overall Convective Heat Transfer Coefficient [J/m^2.K.s]') 
kpa = model.kpa = Var(l, domain=PositiveReals,
    doc='Thermal Conductivity of Bed at Minimum Fluidisation [J/m.K.s]') 
Pr = model.Pr = Var(l, domain=PositiveReals,
    doc='Prandlt Number [-]')
tau = model.tau = Var(l, domain=PositiveReals,
    doc='Average Residence Time of Dense PAckets at HX Surface [s]')

# HX Tube Heat Transfer
Q = model.Q = Var(doc='Total Heat Transfered to Heat Exchanger Fluid [J/s]')
Qhx = model.Qhx = Var(l,
                doc='Heat Transfered to HX Fluid in Discrete Element [J/s]')

# -----------------------------------------------------------------------------
# Outlet Variables

# Gas Outlet Conditions
Gas_Out_F = model.Gas_Out_F = Var(domain=PositiveReals,
                    doc='Gas Flowrate at Gas Outlet [mol/s]')
Gas_Out_P = model.Gas_Out_P = Var(domain=PositiveReals,
                    doc='Pressure at Gas Outlet [Pa]')
Gas_Out_T = model.Gas_Out_T = Var(domain=PositiveReals,
                    doc='Gas Temeprature at Gas Outlet [K]')
Gas_Out_y = model.Gas_Out_y = Var(GasList, domain=PositiveReals,
                    doc='Gas Mole Fractions at Gas Outlet [mol/mol]')

# Solids Inlet Conditions
Solid_Out_F = model.Solid_Out_F = Var(domain=PositiveReals,
                    doc='Solid Flowrate at Solid Outlet [kg/s]')
Solid_Out_T = model.Solid_Out_T = Var(domain=PositiveReals,
                    doc='Solid Temperature at Solid Outlet [K]')
Solid_Out_w = model.Solid_Out_w = Var(SolidList, domain=PositiveReals,
                    doc='Solid Loading at Solid Outlet [mol/kg]')

# HX Fluid Inlet Conditions
HX_Out_F = model.HX_Out_F = Var(domain=PositiveReals,
                    doc='HX Fluid Flowrate at Outlet [mol/s]')
HX_Out_T = model.HX_Out_T = Var(domain=PositiveReals,
                    doc='HX Fluid Temperature at Outlet [K]')
HX_Out_P = model.HX_Out_P = Var(domain=PositiveReals,
                    doc='HX Fluid Pressure at Outlet [Pa]')
HX_Out_y = model.HX_Out_y = Var(HXList, domain=PositiveReals,
                    doc='HX Fluid Mole fractions at Outlet [mol/mol]')

# -----------------------------------------------------------------------------
# Temporary Property Variables

# Gas Phase Properties
V_in = model.V_in = Var(domain=PositiveReals)                     # m^3/mol
V = model.V = Var(l, domain=PositiveReals)                     # m^3/mol
cpg_mol = model.cpgmol = Var(domain = PositiveReals)           # J/mol.K
D = model.D = Var(GasList, l, domain=PositiveReals)            # m^2/s
kg = model.kg = Var(domain = PositiveReals)                    # W/m.K
rhog = model.rhog = Var(l, domain=PositiveReals)               # kg/m^3
mug = model.mug = Var(domain=PositiveReals)                    # Pa.s
MW = model.MW = Var(l, domain=PositiveReals)                   # kg/mol

cpgcsb = model.cpgcsb = Var(GasList, domain = PositiveReals)                    # J/mol.K
cpgcst = model.cpgcst = Var(GasList, domain = PositiveReals)                    # J/mol.K
cpgcsc = model.cpgcsc = Var(GasList, domain = PositiveReals)                    # J/mol.K
cpgcse = model.cpgcse = Var(GasList, domain = PositiveReals)                    # J/mol.K
cpgcgc = model.cpgcgc = Var(GasList, domain = PositiveReals)                    # J/mol.K
cpgcge = model.cpgcge = Var(GasList, domain = PositiveReals)                    # J/mol.K

ap = model.ap = Var()
cps = model.cps = Var()                    # J/kg.K
dp = model.dp = Var()
emf = model.emf = Var()
kp = model.kp = Var()                       # W/m.K
phis = model.phis = Var()
rhos = model.rhos = Var()                   # kg/m^3
F = model.F = Var(domain=NonNegativeReals)

A1 = model.A1 = Var()
A2 = model.A2 = Var()
A3 = model.A3 = Var()
dH1 = model.dH1 = Var()
dH2 = model.dH2 = Var()
dH3 = model.dH3 = Var()
dS1 = model.dS1 = Var()
dS2 = model.dS2 = Var()
dS3 = model.dS3 = Var()
E1 = model.E1 = Var()
E2 = model.E2 = Var()
E3 = model.E3 = Var()
nv = model.nv = Var()

k1c = model.k1c = Var(l, domain=NonNegativeReals)
k2c = model.k2c = Var(l, domain=NonNegativeReals)
k3c = model.k3c = Var(l, domain=NonNegativeReals)
k1e = model.k1e = Var(l, domain=NonNegativeReals)
k2e = model.k2e = Var(l, domain=NonNegativeReals)
k3e = model.k3e = Var(l, domain=NonNegativeReals)

K1c = model.K1c = Var(l, domain=NonNegativeReals)
K2c = model.K2c = Var(l, domain=NonNegativeReals)
K3c = model.K3c = Var(l, domain=NonNegativeReals)
K1e = model.K1e = Var(l, domain=NonNegativeReals)
K2e = model.K2e = Var(l, domain=NonNegativeReals)
K3e = model.K3e = Var(l, domain=NonNegativeReals)

r1c = model.r1c = Var(l, domain=Reals)
r2c = model.r2c = Var(l, domain=Reals)
r3c = model.r3c = Var(l, domain=Reals)
r1e = model.r1e = Var(l, domain=Reals)
r2e = model.r2e = Var(l, domain=Reals)
r3e = model.r3e = Var(l, domain=Reals)

rgc = model.rgc = Var(GasList, l)
rge = model.rge = Var(GasList, l)
rsc = model.rsc = Var(SolidList, l)
rse = model.rse = Var(SolidList, l)

hs_in = model.hs_in = Var(within=Reals)                    # J/kg
hsc = model.hsc = Var(l, within=Reals)                    # J/kg
hse = model.hse = Var(l, within=Reals)                    # J/kg

hhx_in = model.hhx_in = Var()                               # J/mol
hhx = model.hhx = Var(l)                                    # J/mol
rhohx = model.rhohx = Var(domain = PositiveReals)

# =============================================================================
# Properties

model.eq_q1 = Constraint(expr = Gas_In_P*V_in == R*Gas_In_T)

def rule_eq_q2(m, i):
    return P[i]*V[i] == R*Tgb[i]
model.eq_q2 = Constraint(l, rule=rule_eq_q2)

def rule_eq_q3(m, i):
    return MW[i] == (ye['CO2',i]*44 + ye['H2O',i]*18 + ye['N2',i]*28)*1e-3
model.eq_q3 = Constraint(l, rule=rule_eq_q3)

def rule_eq_q4(m, i):
    return rhog[i] == MW[i]*P[i]/(R*Tge[i])
model.eq_q4 = Constraint(l, rule=rule_eq_q4)

def rule_eq_q5(m, i, j):
    if j == 'CO2':
        return D[j,i] == 1e-4*((0.1593-0.1282*(P[i]*1e-5-1.4) + 0.001*(Tge[i]-333.15)+0.0964*((P[i]*1e-5-1.4)**2) -
                           0.0006921*((P[i]*1e-5-1.4)*(Tge[i]-333.15)) -
                           3.3532e-06*(Tge[i]-333.15)**2)*ye['H2O',i]/(ye['H2O',i]+ye['N2',i]) + \
                           (0.1495-0.1204*(P[i]*1e-5-1.4)+0.0008896*(Tge[i]-333.15)+0.0906*((P[i]*1e-5-1.4)**2) -
                            0.0005857*(P[i]*1e-5-1.4)*(Tge[i]-333.15) -
                            3.559e-06*(Tge[i]-333.15)**2)*ye['N2',i]/(ye['H2O',i]+ye['N2',i]))
    elif j == 'H2O':
        return D[j,i] == 1e-4*((0.1593-0.1282*(P[i]*1e-5-1.4)+0.001*(Tge[i]-333.15) +
                           0.0964*((P[i]*1e-5-1.4)**2)-0.0006921*((P[i]*1e-5-1.4)*(Tge[i]-333.15)) -
                           3.3532e-06*(Tge[i]-333.15)**2)*ye['CO2',i]/(ye['CO2',i]+ye['N2',i]) +\
                          (0.2165-0.1743*(P[i]*1e-5-1.4)+0.001377*(Tge[i]-333.15)+0.13109*((P[i]*1e-5-1.4)**2) -
                           0.0009115*(P[i]*1e-5-1.4)*(Tge[i]-333.15) -
                           4.8394e-06*(Tge[i]-333.15)**2)*ye['N2',i]/(ye['CO2',i]+ye['N2',i]))
    else:
        return D[j,i] == 1e-4*((0.1495-0.1204*(P[i]*1e-5-1.4)+0.0008896*(Tge[i]-333.15)+0.0906*((P[i]*1e-5-1.4)**2) -
                           0.0005857*(P[i]*1e-5-1.4)*(Tge[i]-333.15) -
                           3.559e-06*(Tge[i]-333.15)**2)*ye['CO2',i]/(ye['H2O',i]+ye['CO2',i]) +\
                          (0.2165-0.1743*(P[i]*1e-5-1.4)+0.001377*(Tge[i]-333.15)+0.13109*((P[i]*1e-5-1.4)**2) -
                           0.0009115*(P[i]*1e-5-1.4)*(Tge[i]-333.15) -
                           4.8394e-06*(Tge[i]-333.15)**2)*ye['H2O',i]/(ye['H2O',i]+ye['CO2',i]))
model.eq_q5 = Constraint(l, GasList, rule=rule_eq_q5)

cpg_mol.fix(30.0912)
kg.fix(0.0280596)
mug.fix(1.92403e-5)

cpgcsb['CO2'].fix(38.7)
cpgcsb['H2O'].fix(34.1)
cpgcsb['N2'].fix(29.2)

cpgcst['CO2'].fix(38.7)
cpgcst['H2O'].fix(34.1)
cpgcst['N2'].fix(29.2)

cpgcsc['CO2'].fix(38.7)
cpgcsc['H2O'].fix(34.1)
cpgcsc['N2'].fix(29.2)

cpgcse['CO2'].fix(38.7)
cpgcse['H2O'].fix(34.1)
cpgcse['N2'].fix(29.2)

cpgcgc['CO2'].fix(38.7)
cpgcgc['H2O'].fix(34.1)
cpgcgc['N2'].fix(29.2)

cpgcge['CO2'].fix(38.7)
cpgcge['H2O'].fix(34.1)
cpgcge['N2'].fix(29.2)

# Solid Properties
ap.fix(90.49773755656109)
cps.fix(1130)
dp.fix(1.5e-4)
emf.fix(0.5)
kp.fix(1.36)
phis.fix(1)
rhos.fix(442)
F.fix(0)

A1.fix(0.17583)
A2.fix(0.091098)
A3.fix(141.944)
dH1.fix(-72580.3)
dH2.fix(-77079)
dH3.fix(-109691)
dS1.fix(-141.425)
dS2.fix(-216.244)
dS3.fix(-281.255)
E1.fix(29622.8)
E2.fix(83173.5)
E3.fix(27522.5)
nv.fix(1900.46)

def rule_eq_r1(m, i):
    return k1c[i] == A1*Tsc[i]*exp(-E1/(R*Tsc[i]))
model.eq_r1 = Constraint(l, rule=rule_eq_r1)

def rule_eq_r2(m, i):
    return k2c[i] == A2*Tsc[i]*exp(-E2/(R*Tsc[i]))
model.eq_r2 = Constraint(l, rule=rule_eq_r2)

def rule_eq_r3(m, i):
    return k3c[i] == A3*Tsc[i]*exp(-E3/(R*Tsc[i]))
model.eq_r3 = Constraint(l, rule=rule_eq_r3)

def rule_eq_r4(m, i):
    return k1e[i] == A1*Tse[i]*exp(-E1/(R*Tse[i]))
model.eq_r4 = Constraint(l, rule=rule_eq_r4)

def rule_eq_r5(m, i):
    return k2e[i] == A2*Tse[i]*exp(-E2/(R*Tse[i]))
model.eq_r5 = Constraint(l, rule=rule_eq_r5)

def rule_eq_r6(m, i):
    return k3e[i] == A3*Tse[i]*exp(-E3/(R*Tse[i]))
model.eq_r6 = Constraint(l, rule=rule_eq_r6)

def rule_eq_r7(m, i):
    return K1c[i]*P[i] == exp(-dH1/(R*Tsc[i]) + dS1/R)
model.eq_r7 = Constraint(l, rule=rule_eq_r7)

def rule_eq_r8(m, i):
    return K2c[i]*P[i] == exp(-dH2/(R*Tsc[i]) + dS2/R)
model.eq_r8 = Constraint(l, rule=rule_eq_r8)

def rule_eq_r9(m, i):
    return K3c[i]*P[i] == exp(-dH3/(R*Tsc[i]) + dS3/R)
model.eq_r9 = Constraint(l, rule=rule_eq_r9)

def rule_eq_r10(m, i):
    return K1e[i]*P[i] == exp(-dH1/(R*Tse[i]) + dS1/R)
model.eq_r10 = Constraint(l, rule=rule_eq_r10)

def rule_eq_r11(m, i):
    return K2e[i]*P[i] == exp(-dH2/(R*Tse[i]) + dS2/R)
model.eq_r11 = Constraint(l, rule=rule_eq_r11)

def rule_eq_r12(m, i):
    return K3e[i]*P[i] == exp(-dH3/(R*Tse[i]) + dS3/R)
model.eq_r12 = Constraint(l, rule=rule_eq_r12)

def rule_eq_r13(m, i):
    return r1c[i] == k1c[i]*((P[i]*yc['H2O',i]) - (wc['H2O',i]*rhos/K1c[i]))
model.eq_r13 = Constraint(l, rule=rule_eq_r13)

def rule_eq_r14(m, i):
    return r2c[i] == k2c[i]*((1-2*(wc['Car',i]*rhos/nv)-(wc['Bic',i]*rhos/nv))*wc['H2O',i]*rhos*P[i]*yc['CO2',i] -
                                 (((wc['Car',i]*rhos/nv)+(wc['Bic',i]*rhos/nv))*wc['Bic',i]*rhos/K2c[i]))
model.eq_r14 = Constraint(l, rule=rule_eq_r14)

def rule_eq_r15(m, i):
    return r3c[i] == k3c[i]*(((1-2*(wc['Car',i]*rhos/nv)-(wc['Bic',i]*rhos/nv))**2)*(P[i]*yc['CO2',i]) -
                                 ((wc['Car',i]*rhos/nv)*((wc['Car',i]*rhos/nv)+(wc['Bic',i]*rhos/nv))/K3c[i]))
model.eq_r15 = Constraint(l, rule=rule_eq_r15)

def rule_eq_r16(m, i):
    return r1e[i] == k1e[i]*((P[i]*ye['H2O',i]) - (we['H2O',i]*rhos/K1e[i]))
model.eq_r16 = Constraint(l, rule=rule_eq_r16)

def rule_eq_r17(m, i):
    return r2e[i] == k2e[i]*((1-2*(we['Car',i]*rhos/nv)-(we['Bic',i]*rhos/nv))*we['H2O',i]*rhos*P[i]*ye['CO2',i] -
                                (((we['Car',i]*rhos/nv)+(we['Bic',i]*rhos/nv))*we['Bic',i]*rhos/K2e[i]))
model.eq_r17 = Constraint(l, rule=rule_eq_r17)

def rule_eq_r18(m, i):
    return r3e[i] == k3e[i]*(((1-2*(we['Car',i]*rhos/nv)-(we['Bic',i]*rhos/nv))**2)*(P[i]*ye['CO2',i]) -
                                 ((we['Car',i]*rhos/nv)*((we['Car',i]*rhos/nv)+(we['Bic',i]*rhos/nv))/K3e[i]))
model.eq_r18 = Constraint(l, rule=rule_eq_r18)

def rule_eq_r19(m, i, j):
    if j == 'CO2':
        return rgc[j,i] == (nv*r3c[i] + r2c[i])*dl*Ax*delta[i]*fcw[i]*(1-ed[i])
    elif j == 'H2O':
        return rgc[j,i] == r1c[i]*dl*Ax*delta[i]*fcw[i]*(1-ed[i])
    else:
        return rgc[j,i] == 0
model.eq_r19 = Constraint(l, GasList, rule=rule_eq_r19)

def rule_eq_r20(m, i, j):
    if j == 'CO2':
        return rge[j,i] == (nv*r3e[i] + r2e[i])*dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])
    elif j == 'H2O':
        return rge[j,i] == r1e[i]*dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])
    else:
        return rge[j,i] == 0
model.eq_r20 = Constraint(l, GasList, rule=rule_eq_r20)

def rule_eq_r21(m, i, j):
    if j == 'Bic':
        return rsc[j,i] == r2c[i]*dl*Ax*delta[i]*fcw[i]*(1-ed[i])
    elif j == 'H2O':
        return rsc[j,i] == (r1c[i] - r2c[i])*dl*Ax*delta[i]*fcw[i]*(1-ed[i])
    else:
        return rsc[j,i] == nv*r3c[i]*dl*Ax*delta[i]*fcw[i]*(1-ed[i])
model.eq_r21 = Constraint(l, SolidList, rule=rule_eq_r21)

def rule_eq_r22(m, i, j):
    if j == 'Bic':
        return rse[j,i] == r2e[i]*dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])
    elif j == 'H2O':
        return rse[j,i] == (r1e[i] - r2e[i])*dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])
    else:
        return rse[j,i] == nv*r3e[i]*dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])
model.eq_r22 = Constraint(l, SolidList, rule=rule_eq_r22)

model.eq_s1 = Constraint(expr = hs_in == ((Solid_In_w['H2O']+Solid_In_w['Bic'])*(cpgcse['H2O']*Solid_In_T + dH1) +
                        Solid_In_w['Bic']*(cpgcse['CO2']*Solid_In_T + dH2) +
                        Solid_In_w['Car']*(cpgcse['CO2']*Solid_In_T + dH3)) + cps*Solid_In_T)

def rule_eq_s2(m, i):
    return m.hsc[i] == ((wc['H2O',i]+wc['Bic',i])*(cpgcsc['H2O']*Tsc[i] + dH1) +
                        wc['Bic',i]*(cpgcsc['CO2']*Tsc[i] + dH2) +
                        wc['Car',i]*(cpgcsc['CO2']*Tsc[i] + dH3)) + cps*Tsc[i]
model.eq_s2 = Constraint(l, rule=rule_eq_s2)

def rule_eq_s3(m, i):
    return m.hse[i] == ((we['H2O',i]+we['Bic',i])*(cpgcse['H2O']*Tse[i] + dH1) +
                        we['Bic',i]*(cpgcse['CO2']*Tse[i] + dH2) +
                        we['Car',i]*(cpgcse['CO2']*Tse[i] + dH3)) + cps*Tse[i]
model.eq_s3 = Constraint(l, rule=rule_eq_s3)

def rule_eq_t1(m, i):
    return hhx[i] == ((Thx[i]-33.2104)/14170.15 - 0.285)*1e6
model.eq_t1 = Constraint(l, rule=rule_eq_t1)

model.eq_t2 = Constraint(expr = hhx_in == ((HX_In_T-33.2104)/14170.15 - 0.285)*1e6)

rhohx.fix(985.9393497324021)

# =============================================================================
# Model Equations

# -----------------------------------------------------------------------------
# Reactor Sizing and Design

# a1 - Differential Length
model.eq_a1 = Constraint(expr = dl*ND == Lb)

# a2 - Total Reactor Cross-Sectional Area (incl. HX tubes)
model.eq_a2 = Constraint(expr = Areact == Ax + (pi/4)*dx**2*Nx)

# a3, a4 - HX Tube Pitch and Spacing
model.eq_a3 = Constraint(expr = Areact == Nx*lp**2)
model.eq_a4 = Constraint(expr = lhx*10 == lp*10 - dx*10)

# a5 - Reactor Cross-Sectional Area
model.eq_a5 = Constraint(expr = Areact == 0.25*pi*Dt**2)

# a6 - Distributor Design
model.eq_a6 = Constraint(expr = 10 == nor*Ao*10)

# a7 - Surface Area of HX Tubes
model.eq_a7 = Constraint(expr = Ahx == pi*dx*Lb*Nx)

# a8 - Hydraulic Diameter
model.eq_a8 = Constraint(expr = 0.25*Dte*pi*(Dt+dx*Nx) == Ax)

# -----------------------------------------------------------------------------
# Mass and Energy Balances

# b1 - Bubble Region Gas Component Balances
def rule_eq_b1(m, i, j):
    if i == 1:
        return 0 == Gas_In_F*Gas_In_y[j] - Gb[i]*yb[j,i] - dl*Ax*delta[i]*Kbc[j,i]*(cb[j,i]-cc[j,i]) + Kgbulk_c[j,i]
    else:
        return 0 == Gb[i-1]*yb[j,i-1] - Gb[i]*yb[j,i] - dl*Ax*delta[i]*Kbc[j,i]*(cb[j,i]-cc[j,i]) + Kgbulk_c[j,i]
model.eq_b1 = Constraint(l, GasList, rule=rule_eq_b1)

# b2 - Cloud-Wake Region Gas Component Balances
def rule_eq_b2(m, i, j):
    return 0 == dl*Ax*delta[i]*Kbc[j,i]*(cb[j,i] - cc[j,i]) - dl*Ax*delta[i]*Kce[j,i]*(cc[j,i] - ce[j,i]) - rgc[j,i]
model.eq_b2 = Constraint(l, GasList, rule=rule_eq_b2)

# b3 - Emulsion Region Gas Component Balances
def rule_eq_b3(m, i, j):
    return 0 == dl*Ax*delta[i]*Kce[j,i]*(cc[j,i] - ce[j,i]) - Kgbulk_c[j,i] - rge[j,i]
model.eq_b3 = Constraint(l, GasList, rule=rule_eq_b3)

# b4 - Bubble Region Gas Energy Balance
def rule_eq_b4(m, i):
    if i == 1:
        return 0 == Gas_In_F*Gas_In_T*cpg_mol - Ghb[i] - dl*Ax*delta[i]*Hbc[i]*(Tgb[i]-Tgc[i]) + Hgbulk[i]
    else:
        return 0 == Ghb[i-1] - Ghb[i] - dl*Ax*delta[i]*Hbc[i]*(Tgb[i]-Tgc[i]) + Hgbulk[i]
model.eq_b4 = Constraint(l, rule=rule_eq_b4)

# b5 - Cloud-Wake Region Gas Energy Balance
def rule_eq_b5(m, i):
    return 0 == dl*Ax*delta[i]*Hbc[i]*(Tgb[i]-Tgc[i]) - dl*Ax*delta[i]*Hce[i]*(Tgc[i]-Tge[i]) - \
                dl*Ax*delta[i]*fcw[i]*(1-ed[i])*rhos*ap*hp[i]*(Tgc[i]-Tsc[i]) - \
                sum(rgc[j,i]*cpgcgc[j] for j in GasList)*Tgc[i]
model.eq_b5 = Constraint(l, rule=rule_eq_b5)

# b6 - Emulsion Region Gas Energy Balance
def rule_eq_b6(m, i):
    return 0 == dl*Ax*delta[i]*Hce[i]*(Tgc[i]-Tge[i]) - \
                dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])*rhos*ap*hp[i]*(Tge[i]-Tse[i]) -\
                Hgbulk[i] -\
                sum(rge[j,i]*cpgcge[j] for j in GasList)*Tge[i]
model.eq_b6 = Constraint(l, rule=rule_eq_b6)

# b7 - Cloud-Wake Region Sorbent Balance
def rule_eq_b7(m, i):
    if i == 1:
        return 0 == -Jc[i]*Ax + Ksbulk[i]
    else:
        return 0 == Jc[i-1]*Ax - Jc[i]*Ax + Ksbulk[i]
model.eq_b7 = Constraint(l, rule=rule_eq_b7)

# b8 - Emulsion Region Sorbent Balance
def rule_eq_b8(m, i):
    if (i == ND) or (i == 1):
        return Constraint.Skip
    else:
        return 0 == Je[i+1]*Ax - Je[i]*Ax - Ksbulk[i]
model.eq_b8 = Constraint(l, rule=rule_eq_b8)

# b9 - Cloud-Wake Region Adsorbed Species Balance
def rule_eq_b9(m, i, j):
    if i == 1:
        return 0 == Ksbulk_c[j,i] - Jc[i]*Ax*wc[j,i] + rsc[j,i] - Kcebs_c[j,i]*(wc[j,i]-we[j,i])
    else:
        return 0 == Jc[i-1]*Ax*wc[j,i-1] - Jc[i]*Ax*wc[j,i] + Ksbulk_c[j,i] + rsc[j,i] - Kcebs_c[j,i]
model.eq_b9 = Constraint(l, SolidList, rule=rule_eq_b9)

# b10 - Emulsion Region Adsorbed Species Balance
def rule_eq_b10(m, i, j):
    if (i == ND) or (i == 1):
        return Constraint.Skip
    else:
        return 0 == Je[i+1]*Ax*we[j,i+1] - Je[i]*Ax*we[j,i] - Ksbulk_c[j,i] + rse[j,i] + Kcebs_c[j,i]
model.eq_b10 = Constraint(l, SolidList, rule=rule_eq_b10)

# b11 - Cloud-Wake Region Solid Energy Balance
def rule_eq_b11(m, i):
    if i == 1:
        return 0 == -Jhc[1]*Ax + Hsbulk[1] -\
                    dl*Ax*delta[1]*rhos*Kcebs[1]*(hsc[1] - hse[1]) +\
                    dl*Ax*delta[1]*fcw[1]*(1-ed[1])*rhos*ap*hp[1]*(Tgc[1]-Tsc[1]) +\
                    sum((rgc[j,1]*cpgcgc[j]) for j in GasList)*(Tgc[1])
    else:
        return 0 == Jhc[i-1]*Ax - Jhc[i]*Ax + Hsbulk[i] -\
                    dl*Ax*delta[i]*rhos*Kcebs[i]*(hsc[i] - hse[i]) +\
                    dl*Ax*delta[i]*fcw[i]*(1-ed[i])*rhos*ap*hp[i]*(Tgc[i]-Tsc[i]) +\
                    sum((rgc[j,i]*cpgcgc[j]) for j in GasList)*(Tgc[i])
model.eq_b11 = Constraint(l, rule=rule_eq_b11)

# b12 - Emulsion Region Solid Energy Balance
def rule_eq_b12(m, i):
    if (i == ND) or (i == 1):
        return Constraint.Skip
    else:
        return 0 == Jhe[i+1]*Ax - Jhe[i]*Ax - Hsbulk[i] +\
                    dl*Ax*delta[i]*rhos*Kcebs[i]*(hsc[i]-hse[i]) +\
                    dl*Ax*(1-fcw[i]*delta[i]-delta[i])*(1-ed[i])*rhos*ap*hp[i]*(Tge[i]-Tse[i]) +\
                    sum((rge[j,i]*cpgcge[j]) for j in GasList)*Tge[i] +\
                    Qhx[i]
model.eq_b12 = Constraint(l, rule=rule_eq_b12)

# -----------------------------------------------------------------------------
# Solid Phase Boundary Conditions (Emulsion Region)

# Bottom Feed, Underflow Outlet
model.eq_b13 = Constraint(expr = 0 == Je[2]*Ax + Solid_In_F - Ksbulk[1] - Je[1]*Ax)
model.eq_b14 = Constraint(expr = 0 == Jc[value(ND)]*Ax - Je[value(ND)]*Ax - Ksbulk[value(ND)])

def rule_eq_b15(m, j):
    return 0 == Je[2]*Ax*we[j,2] + Solid_In_F*Solid_In_w[j] - Je[1]*Ax*we[j,1] - Ksbulk_c[j,1] + rse[j,1] + Kcebs_c[j,1]
model.eq_b15 = Constraint(SolidList, rule=rule_eq_b15)

def rule_eq_b16(m, j):
    return 0 == Jc[value(ND)]*Ax*wc[j,value(ND)] - Je[value(ND)]*Ax*we[j,value(ND)] - Ksbulk_c[j,value(ND)] + rse[j,value(ND)] + Kcebs_c[j,value(ND)]
model.eq_b16 = Constraint(SolidList, rule=rule_eq_b16)

model.eq_b17 = Constraint(expr = 0 == Jhe[2]*Ax + Solid_In_F*hs_in - Jhe[1]*Ax - Hsbulk[1] +\
                    dl*Ax*delta[1]*rhos*Kcebs[1]*(hsc[1]-hse[1]) +\
                    dl*Ax*(1-fcw[1]*delta[1]-delta[1])*(1-ed[1])*rhos*ap*hp[1]*(Tge[1]-Tse[1]) +\
                    sum((rge[j,1]*cpgcge[j]) for j in GasList)*Tge[1] +\
                    Qhx[1])

model.eq_b18 = Constraint(expr = 0 == Jhc[value(ND)]*Ax - Jhe[value(ND)]*Ax - Hsbulk[value(ND)] +\
                    dl*Ax*delta[value(ND)]*rhos*Kcebs[value(ND)]*(hsc[value(ND)]-hse[value(ND)]) +\
                    dl*Ax*(1-fcw[value(ND)]*delta[value(ND)]-delta[value(ND)])*(1-ed[value(ND)])*rhos*ap*hp[value(ND)]*(Tge[value(ND)]-Tse[value(ND)]) +\
                    sum((rge[j,value(ND)]*cpgcge[j]) for j in GasList)*Tge[value(ND)] +\
                    Qhx[value(ND)])

# -----------------------------------------------------------------------------
# Flowrate and Flux Relationships

# c1 - Superficial Gas Velocity
def rule_eq_c1(m, i):
    return vg[i]*Ax == Gb[i]*V[i]
model.eq_c1 = Constraint(l, rule=rule_eq_c1)

# c2 - Bubble Gas Flowrate
def rule_eq_c2(m, i): #a13
    return vg[i] == vb[i]*delta[i]
model.eq_c2 = Constraint(l, rule=rule_eq_c2)

# c3 - Bubble Gas Enthalpy Flowrate
def rule_eq_c3(m, i):
    return Ghb[i] == Gb[i]*Tgb[i]*cpg_mol
model.eq_c3 = Constraint(l, rule=rule_eq_c3)

# c4 - Cloud-Wake Region Solids Flux
def rule_eq_c4(m, i):
    return Jc[i] == fw*delta[i]*rhos*(1-ed[i])*vb[i]
model.eq_c4 = Constraint(l, rule=rule_eq_c4)

# c5 - CW Solids Enthalpy Flux
def rule_eq_c5(m, i):
    return Jhc[i] == Jc[i]*hsc[i]
model.eq_c5 = Constraint(l, rule=rule_eq_c5)

# c6 - Emulsion Region Solids Velocity
def rule_eq_c6(m, i):
    return Je[i] == (1-fcw[i]*delta[i]-delta[i])*rhos*(1-ed[i])*us[i]
model.eq_c6 = Constraint(l, rule=rule_eq_c6)

# c7 - CW Solids Enthalpy Flux
def rule_eq_c7(m, i):
    return Jhe[i] == Je[i]*hse[i]
model.eq_c7 = Constraint(l, rule=rule_eq_c7)

# -----------------------------------------------------------------------------
# Material Inlets and Outlets

# d1-d4 -  Gas Outlet
model.eq_d1 = Constraint(expr = Gas_Out_F == Gb[value(ND)])
model.eq_d2 = Constraint(expr = Gas_Out_P == P[value(ND)])
model.eq_d3 = Constraint(expr = Gas_Out_T == Tgb[value(ND)])
def rule_eq_d4(m, j):
    return Gas_Out_y[j] == yb[j,value(ND)]
model.eq_d4 = Constraint(GasList, rule=rule_eq_d4)

# d5-d7 - Solid Outlets                                                             Change for different boundary conditions
model.eq_d5 = Constraint(expr = Solid_Out_F == Je[1]*Ax)
model.eq_d6 = Constraint(expr = Solid_Out_T == Tse[1])
def rule_eq_d7(m, j):
    return Solid_Out_w[j] == we[j,1]
model.eq_d7 = Constraint(SolidList, rule=rule_eq_d7)

# d8-d11 - HX Fluid Outlet Conditions
model.eq_d8 = Constraint(expr = HX_Out_F == HX_In_F)
model.eq_d9 = Constraint(expr = HX_Out_P == Phx[1])
model.eq_d10 = Constraint(expr = HX_Out_T == Thx[1])
def rule_eq_d11(m, j):
    return HX_Out_y[j] == HX_In_y[j]
model.eq_d11 = Constraint(HXList, rule=rule_eq_d11)

# -----------------------------------------------------------------------------
# Mole Fraction Relationships

def rule_eq_e1(m, i, j):
    return cb[j, i]*V[i] == yb[j, i]
model.eq_e1 = Constraint(l, GasList, rule=rule_eq_e1)

def rule_eq_e2(m, i, j):
    return cc[j, i] == yc[j, i]*cct[i]
model.eq_e2 = Constraint(l, GasList, rule=rule_eq_e2)

def rule_eq_e3(m, i, j):
    return ce[j, i] == ye[j, i]*cet[i]
model.eq_e3 = Constraint(l, GasList, rule=rule_eq_e3)

def rule_eq_e4(m, i):
    return 1 == sum(cb[j, i] for j in GasList)*V[i]
model.eq_e4 = Constraint(l, rule=rule_eq_e4)

def rule_eq_e5(m, i):
    return cct[i] == sum(cc[j, i] for j in GasList)
model.eq_e5 = Constraint(l, rule=rule_eq_e5)

def rule_eq_e6(m, i):
    return cet[i] == sum(ce[j, i] for j in GasList)
model.eq_e6 = Constraint(l, rule=rule_eq_e6)

# -----------------------------------------------------------------------------
# Bulk Flow and Mixing Relationships

# f1 - Bulk Gas Mass Transfer
def rule_eq_f1(m, i, j):
    return Kgbulk_c[j,i] == (6*Kd*delta[i]*dl*Ax*(cet[i]-(1/V[i]))/db[i])*yb[j,i]    # Semi-updated version
model.eq_f1 = Constraint(l, GasList, rule=rule_eq_f1)

# f2 - Bulk Gas Energy Transfer
def rule_eq_f2(m, i):
    return Hgbulk[i] == (6*Kd*delta[i]*dl*Ax*(cet[i]-(1/V[i]))/db[i])*cpg_mol*Tgb[i]    # Semi-updated version
model.eq_f2 = Constraint(l, rule=rule_eq_f2)

# f3 - Bulk Solids Mass Transfer
def rule_eq_f3(m, i, j):
    return Ksbulk_c[j,i] == Ksbulk[i]*we[j,i]
model.eq_f3 = Constraint(l, SolidList, rule=rule_eq_f3)

# f4 - Bulk Solids Energy Transfer
def rule_eq_f4(m, i):
    return Hsbulk[i] == Ksbulk[i]*hse[i]
model.eq_f4 = Constraint(l, rule=rule_eq_f4)

# f5 - Bulk Solids Component Mixing
def rule_eq_f5(m, i, j):
    return Kcebs_c[j,i] == dl*Ax*delta[i]*rhos*Kcebs[i]*(wc[j,i]-we[j,i])
model.eq_f5 = Constraint(l, SolidList, rule=rule_eq_f5)

# -----------------------------------------------------------------------------
# g1 - Pressure Drop
def rule_eq_g1(m, i):
    if i == 1:
        return P[i] - Gas_In_P + 3400 == -dl*(1-e[i])*rhos*gc
    else:
        return P[i] - P[i-1] == -dl*(1-e[i])*rhos*gc
model.eq_g1 = Constraint(l, rule=rule_eq_g1)

# -----------------------------------------------------------------------------
# Fluidisation Conditions and Bubble Behaviour

# h1 - Archimedes Number
def rule_eq_h1(m, i):
    return Ar[i]*mug**2 == (dp**3)*rhog[i]*(rhos-rhog[i])*gc
model.eq_h1 = Constraint(l, rule=rule_eq_h1)

# h2 - Emulsion Region Gas Velocity
def rule_eq_h2(m, i):
    return ve[i]*((dp**0.568)*(gc**0.663)*((rhos-rhog[i])**0.663)*((i*dl)**0.244)) == vmf*188*(rhog[i]**0.089)*(mug**0.371)*exp(0.508*F);
model.eq_h2 = Constraint(l, rule=rule_eq_h2)

# h3 - Initial Bubble Diameter
model.eq_h3 = Constraint(expr = db0 == 1.38*(gc**(-0.2))*(((Gas_In_F*V_in/Ax)-ve[1])*Ao)**0.4)

# h4 - Bubble Size Coefficients
model.eq_h4 = Constraint(expr = (g1*vmf)**2 == (2.56E-2**2)*(Dt/gc))

# h5 - Maximum Bubble Diameter
def rule_eq_h5(m, i):
    return (dbm[i]**5)*gc == (2.59**5)*((vg[i]-ve[i])*Ax)**2
model.eq_h5 = Constraint(l, rule=rule_eq_h5)

# h6 - Constrained Bubble Diameter
def rule_eq_h6(m, i):
    if i == 1:
        return (db[i] - db0) == 0.3*(dl/Dt)*(dbm[i] - db[i] - g1*(Dt*db[i])**0.5)
    else:
        return (db[i] - db[i-1]) == 0.3*(dl/Dt)*(dbm[i] - db[i] - g1*(Dt*db[i])**0.5)
model.eq_h6 = Constraint(l, rule=rule_eq_h6)

# h7 - Bubble Rise Velocity
def rule_eq_h7(m, i):
    return vbr[i]**2 == (0.711**2)*(gc*db[i])
model.eq_h7 = Constraint(l, rule=rule_eq_h7)

# h8 - Bubble Velocity
def rule_eq_h8(m, i):
    return vb[i] == 1.55*((vg[i]-vmf)+14.1*(db[i]+0.005))*(Dte**0.32) + vbr[i]
model.eq_h8 = Constraint(l, rule=rule_eq_h8)

# h9 - Cloud to Bubble Volume Ratio
def rule_eq_h9(m, i):
    return fc[i]*(vbr[i]-(vmf/emf)) == 3*(vmf/emf)
model.eq_h9 = Constraint(l, rule=rule_eq_h9)

# h10 - Cloud-Wake to Bubble Volume Ratio
def rule_eq_h10(m, i):
    return fcw[i] == fc[i] + fw
model.eq_h10 = Constraint(l, rule=rule_eq_h10)

# h11 - Emulsion Region Voidage
def rule_eq_h11(m, i):
    return (1 - emf)*((dp**0.1)*(gc**0.118)*((rhos-rhog[i])**0.118)*((i*dl)**0.043)) == 2.54*(rhog[i]**0.016)*(mug**0.066)*exp(0.090*F)*(1-ed[i])
model.eq_h11 = Constraint(l, rule=rule_eq_h11)

# h12 - Average Voidage
def rule_eq_h12(m, i):
    return (1 - e[i]) == (1 - ed[i])*(1-delta[i])
model.eq_h12 = Constraint(l, rule=rule_eq_h12)

# -----------------------------------------------------------------------------
# K-L Model Heat and Mass Transfer Coefficients

# i1 - Bubble to Cloud-Wake Gas Mass Transfer Coefficient
def rule_eq_i1(m, i, j):
    return (Kbc[j,i]**4)*(db[i]**5) == ((1.32*4.5*vmf)**4)*db[i] + (5.85**4)*(D[j,i]**2)*gc
model.eq_i1 = Constraint(l, GasList, rule=rule_eq_i1)

# i2 - Cloud-Wake to Emuslion Gas Mass Transfer Coefficient
def rule_eq_i2(m, i, j):
    return Kce[j,i]**2 == (6.77**2)*ed[i]*D[j,i]*vbr[i]/(db[i]**3)
model.eq_i2 = Constraint(l, GasList, rule=rule_eq_i2)

# i3 - Cloud-Wake to Emulsion Solid Mass Transfer Coefficient
def rule_eq_i3(m, i):
    return Kcebs[i]*((1-delta[i])*ed[i]*db[i]) == 3*(1-ed[i])*ve[i]
model.eq_i3 = Constraint(l, rule=rule_eq_i3)

# i4 - Bubble to Cloud-Wake Gas Heat Transfer Coefficient
def rule_eq_i4(m, i):
    return Hbc[i]*db[i]**(5/4) == 1.32*4.5*vmf*cpg_mol*db[i]**(1/4)/V[i] +\
                       5.85*(kg*cpg_mol/V[i])**0.5*(gc**0.25)
model.eq_i4 = Constraint(l, rule=rule_eq_i4)

# i5 - Cloud-Wake to Emulsion Gas Heat Transfer Coefficient
def rule_eq_i5(m, i): #a42
    return Hce[i]**2*db[i]**3 == 6.78**2*ed[i]*vbr[i]*kg*cct[i]*cpg_mol
model.eq_i5 = Constraint(l, rule=rule_eq_i5)

# i6 - Convective Heat Transfer Coefficient
def rule_eq_i6(m, i):
    return 0.03*kg*(ve[i]*dp*rhog[i]/mug)**1.3 == hp[i]*dp
model.eq_i6 = Constraint(l, rule=rule_eq_i6)

# -----------------------------------------------------------------------------
# Mickley and Fairbanks HX Model

# j1 - Thermal Conductivity of Bed at Minimum Fluidisation
def rule_eq_j1(m, i):
    return kpa[i] == (3.58-2.5*ed[i])*kg*((kp/kg)**(0.46-0.46*ed[i]))
model.eq_j1 = Constraint(l, rule=rule_eq_j1)

# j2 - Fluidisation Number
def rule_eq_j2(m, i):
    return fn[i]*vmf == vg[i]
model.eq_j2 = Constraint(l, rule=rule_eq_j2)

# j3 - Residence Time of Emulsion Packets at HX Surface
def rule_eq_j3(m, i):
    return tau[i] == 0.44*((dp*gc/((vmf**2)*((fn[i]-ah)**2)))**0.14)*((dp/dx)**0.225)
model.eq_j3 = Constraint(l, rule=rule_eq_j3)

# j4 - Fraction of Time HX Surface is Exposed to Emulsion Packets
def rule_eq_j4(m, i):
    return fb[i] == 0.33*(((vmf**2)*((fn[i]-ah)**2)/(dp*gc))**0.14)
model.eq_j4 = Constraint(l, rule=rule_eq_j4)

# j5 - Dense Region Heat Transfer Coefficient
def rule_eq_j5(m, i):
    return hd[i]*sqrt(pi*tau[i]) == 2*sqrt(kpa[i]*rhos*cps*(1-ed[i]))
model.eq_j5 = Constraint(l, rule=rule_eq_j5)

# j6 - Bubble Region Heat Transfer Coefficient - Prandlt Number
def rule_eq_j6(m, i):
    return Pr[i]*kg == cpg_mol*mug/MW[i]
model.eq_j6 = Constraint(l, rule=rule_eq_j6)

# j7 - Bubble Region Heat Transfer Coefficient
def rule_eq_j7(m, i):
    return hl[i]*dp == 0.009*(Ar[i]**0.5)*(Pr[i]**0.33)*kg
model.eq_j7 = Constraint(l, rule=rule_eq_j7)

# j8 - Total HX Heat Transfer Coefficient
def rule_eq_j8(m, i):
    return ht[i] == fb[i]*hd[i] + (1-fb[i])*hl[i]
model.eq_j8 = Constraint(l, rule=rule_eq_j8)

# -----------------------------------------------------------------------------
# HX Tube Mass and Energy Balances - Liquid Water Only

# k1 - HX Fluid Pressure Drop
def rule_eq_k1(m, i):
    if i == ND:
        return (Phx[i] - HX_In_P) == rhohx*gc*dl
    else:
        return (Phx[i] - Phx[i+1]) == rhohx*gc*dl
model.eq_k1 = Constraint(l, rule=rule_eq_k1)

# k2 - HX Tube Heat Transfer
def rule_eq_k2(m, i):
    return Qhx[i] == dl*pi*dx*ht[i]*dThx[i]*Nx*Cr
model.eq_k2 = Constraint(l, rule=rule_eq_k2)

# k3 - HX Fluid Energy Balance
def rule_eq_k3(m, i):
    if i == ND:
        return 0 == HX_In_F*(hhx[i]-hhx_in) + Qhx[i]
    else:
        return 0 == HX_In_F*(hhx[i]-hhx[i+1]) + Qhx[i]
model.eq_k3 = Constraint(l, rule=rule_eq_k3)

# k4 - Temperature Difference between Tube and Bed
def rule_eq_k4(m, i):
    return dThx[i] == Ttube[i] - Tse[i]
model.eq_k4 = Constraint(l, rule=rule_eq_k4)

# k5 - HX Tube Wall Energy Balance
def rule_eq_k5(m, i):
    return ht[i]*dThx[i]*Cr == hw*(Thx[i]-Ttube[i])
model.eq_k5 = Constraint(l, rule=rule_eq_k5)

# k6 - Total Heat Duty
model.eq_k6 = Constraint(expr = Q == HX_In_F*(hhx_in-hhx[1]))

# -----------------------------------------------------------------------------
# Particle Model

# Minimum Fluidisation Velocity
model.eq_l1 = Constraint(expr = 1.75/(phis*emf**3)*(dp*vmf*rhog[1]/mug)**2 +
           150*(1-emf)/(phis**2*emf**3)*(dp*vmf*rhog[1]/mug) ==
           dp**3*rhog[1]*(rhos-rhog[1])*gc/mug**2)

# =============================================================================
# Create Connectors

'''model.Gas_In = Connector(model.Gas_In_F, Gas_In_T, Gas_In_P, Gas_In_y)
                    
model.HX_In = Connector(HX_In_F, HX_In_T, HX_In_P, HX_In_y)
                    
model.Solid_In = Connector(Solid_In_F, Solid_In_T, Solid_In_w)
                    
model.Gas_Out = Connector(Gas_Out_F, Gas_Out_T, Gas_Out_P, Gas_Out_y)
                    
model.HX_Out = Connector(HX_Out_F, HX_Out_T, HX_Out_P, HX_Out_y)

model.Solid_Out = Connector(Solid_Out_F, Solid_Out_T, Solid_Out_w)'''

# =============================================================================
# Initialisation Constraints

def rule_IP_b5(m, i):
    return Tgc[i] == Tgb[i]
model.IP_b5 = Constraint(l, rule=rule_IP_b5)

def rule_IP_b6(m, i):
    return Tge[i] == Tgb[i]
model.IP_b6 = Constraint(l, rule=rule_IP_b6)

def rule_IP_b9(m, i, j):
    if j == 'H2O':
        return K1c[i]*P[i]*yc['H2O',i] == wc['H2O',i]*rhos
    elif j == 'Bic':
        return (1-2*(wc['Car',i]*rhos/nv) - (wc['Bic',i]*rhos/nv))*wc['H2O',i]*rhos*P[i]*yc['CO2',i] == (((wc['Car',i]*rhos/nv)+(wc['Bic',i]*rhos/nv))*wc['Bic',i]*rhos/K2c[i])
    else:
        return ((1-2*(wc['Car',i]*rhos/nv)-(wc['Bic',i]*rhos/nv))**2)*(P[i]*yc['CO2',i]) == ((wc['Car',i]*rhos/nv)*((wc['Car',i]*rhos/nv)+(wc['Bic',i]*rhos/nv))/K3c[i])
model.IP_b9 = Constraint(l, SolidList, rule=rule_IP_b9)

def rule_IP_b10(m, i, j):
    if (i == 1) or (i == ND):
        return Constraint.Skip
    else:
        if j == 'H2O':
            return P[i]*ye['H2O',i] == we['H2O',i]*rhos/K1e[i]
        elif j == 'Bic':
            return (1-2*(we['Car',i]*rhos/nv) - (we['Bic',i]*rhos/nv))*we['H2O',i]*rhos*P[i]*ye['CO2',i] == (((we['Car',i]*rhos/nv)+(we['Bic',i]*rhos/nv))*we['Bic',i]*rhos/K2e[i])
        else:
            return ((1-2*(we['Car',i]*rhos/nv)-(we['Bic',i]*rhos/nv))**2)*(P[i]*ye['CO2',i]) == ((we['Car',i]*rhos/nv)*((we['Car',i]*rhos/nv)+(we['Bic',i]*rhos/nv))/K3e[i])
model.IP_b10 = Constraint(l, SolidList, rule=rule_IP_b10)

def rule_IP_b15(m, j):
    if j == 'H2O':
        return P[1]*ye['H2O',1] == we['H2O',1]*rhos/K1e[1]
    elif j == 'Bic':
        return (1-2*(we['Car',1]*rhos/nv) - (we['Bic',1]*rhos/nv))*we['H2O',1]*rhos*P[1]*ye['CO2',1] == (((we['Car',1]*rhos/nv)+(we['Bic',1]*rhos/nv))*we['Bic',1]*rhos/K2e[1])
    else:
        return ((1-2*(we['Car',1]*rhos/nv)-(we['Bic',1]*rhos/nv))**2)*(P[1]*ye['CO2',1]) == ((we['Car',1]*rhos/nv)*((we['Car',1]*rhos/nv)+(we['Bic',1]*rhos/nv))/K3e[1])
model.IP_b15 = Constraint(SolidList, rule=rule_IP_b15)

def rule_IP_b16(m, j):
    if j == 'H2O':
        return P[value(ND)]*ye['H2O',value(ND)] == we['H2O',value(ND)]*rhos/K1e[value(ND)]
    elif j == 'Bic':
        return (1-2*(we['Car',value(ND)]*rhos/nv) - (we['Bic',value(ND)]*rhos/nv))*we['H2O',value(ND)]*rhos*P[value(ND)]*ye['CO2',value(ND)] == (((we['Car',value(ND)]*rhos/nv)+(we['Bic',value(ND)]*rhos/nv))*we['Bic',value(ND)]*rhos/K2e[value(ND)])
    else:
        return ((1-2*(we['Car',value(ND)]*rhos/nv)-(we['Bic',value(ND)]*rhos/nv))**2)*(P[value(ND)]*ye['CO2',value(ND)]) == ((we['Car',value(ND)]*rhos/nv)*((we['Car',value(ND)]*rhos/nv)+(we['Bic',value(ND)]*rhos/nv))/K3e[value(ND)])
model.IP_b16 = Constraint(SolidList, rule=rule_IP_b16)

# =============================================================================
# Fix Value of Controlled Variable

Dt.fix(9)
Lb.fix(5)
nor.fix(2500)
dx.fix(0.02)
lhx.fix(0.24)

Gas_In_F.fix(1613.95)   # mol/s
Gas_In_P.fix(132285)   # Pa
Gas_In_T.fix(316.15)        # K
Gas_In_y['CO2'].fix(0.13)
Gas_In_y['H2O'].fix(0.06)
Gas_In_y['N2'].fix(0.81)

Solid_In_F.fix(165.416) # kg/s
Solid_In_T.fix(338.15)      # K
Solid_In_w['Bic'].fix(0.01)
Solid_In_w['Car'].fix(0.7)
Solid_In_w['H2O'].fix(0.7)

HX_In_F.fix(16666.7)    # mol/s
HX_In_T.fix(306.15)         # s
HX_In_P.fix(112000)     # Pa
HX_In_y['H2O'].fix(1)

# =============================================================================
# Set Initial Values of State and Property Variables

for i in l:
    model.Gb[i] = value(Gas_In_F)
    model.P[i] = value(Gas_In_P)
    model.Phx[i] = value(HX_In_P)
    model.Tgb[i] = value(Gas_In_T)
    model.Tgc[i] = value(Gas_In_T)
    model.Tge[i] = value(Gas_In_T)
    model.Tsc[i] = value(Solid_In_T)
    model.Tse[i] = value(Solid_In_T)
    model.Thx[i] = value(HX_In_T)

    for j in GasList:
        model.yb[j,i] = value(Gas_In_y[j])
        model.yc[j,i] = value(Gas_In_y[j])
        model.ye[j,i] = value(Gas_In_y[j])

# Property Initial Values
for i in l:
    model.V[i] = value(R)*value(Tgb[i])/value(P[i])
    model.MW[i] = (value(ye['CO2',i])*44 + value(ye['H2O',i])*18 + value(ye['N2',i])*28)*1e-3
    model.rhog[i] = value(MW[i])/value(V[i])
    
    for j in GasList:
        if j == 'CO2':
            model.D[j,i] = 1e-4*((0.1593-0.1282*(value(P[i])*1e-5-1.4) + 0.001*(value(Tge[i])-333.15)+0.0964*((value(P[i])*1e-5-1.4)**2) -
                           0.0006921*((value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15)) -
                           3.3532e-06*(value(Tge[i])-333.15)**2)*value(ye['H2O',i])/(value(ye['H2O',i])+value(ye['N2',i])) + \
                           (0.1495-0.1204*(value(P[i])*1e-5-1.4)+0.0008896*(value(Tge[i])-333.15)+0.0906*((value(P[i])*1e-5-1.4)**2) -
                            0.0005857*(value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15) -
                            3.559e-06*(value(Tge[i])-333.15)**2)*value(ye['N2',i])/(value(ye['H2O',i])+value(ye['N2',i])))
        elif j == 'H2O':
            model.D[j,i] = 1e-4*((0.1593-0.1282*(value(P[i])*1e-5-1.4)+0.001*(value(Tge[i])-333.15) +
                           0.0964*((value(P[i])*1e-5-1.4)**2)-0.0006921*((value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15)) -
                           3.3532e-06*(value(Tge[i])-333.15)**2)*value(ye['CO2',i])/(value(ye['CO2',i])+value(ye['N2',i])) +\
                          (0.2165-0.1743*(value(P[i])*1e-5-1.4)+0.001377*(value(Tge[i])-333.15)+0.13109*((value(P[i])*1e-5-1.4)**2) -
                           0.0009115*(value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15) -
                           4.8394e-06*(value(Tge[i])-333.15)**2)*value(ye['N2',i])/(value(ye['CO2',i])+value(ye['N2',i])))
        else:
            model.D[j,i] = 1e-4*((0.1495-0.1204*(value(P[i])*1e-5-1.4)+0.0008896*(value(Tge[i])-333.15)+0.0906*((value(P[i])*1e-5-1.4)**2) -
                           0.0005857*(value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15) -
                           3.559e-06*(value(Tge[i])-333.15)**2)*value(ye['CO2',i])/(value(ye['H2O',i])+value(ye['CO2',i])) +\
                          (0.2165-0.1743*(value(P[i])*1e-5-1.4)+0.001377*(value(Tge[i])-333.15)+0.13109*((value(P[i])*1e-5-1.4)**2) -
                           0.0009115*(value(P[i])*1e-5-1.4)*(value(Tge[i])-333.15) -
                           4.8394e-06*(value(Tge[i])-333.15)**2)*value(ye['H2O',i])/(value(ye['H2O',i])+value(ye['CO2',i])))

    model.k1c[i] = value(A1)*value(Tsc[i])*exp(-value(E1)/(value(R)*value(Tsc[i])))
    model.k2c[i] = value(A2)*value(Tsc[i])*exp(-value(E2)/(value(R)*value(Tsc[i])))
    model.k3c[i] = value(A3)*value(Tsc[i])*exp(-value(E3)/(value(R)*value(Tsc[i])))
    model.k1e[i] = value(A1)*value(Tse[i])*exp(-value(E1)/(value(R)*value(Tse[i])))
    model.k2e[i] = value(A2)*value(Tse[i])*exp(-value(E2)/(value(R)*value(Tse[i])))
    model.k3e[i] = value(A3)*value(Tse[i])*exp(-value(E3)/(value(R)*value(Tse[i])))

    model.K1c[i] == exp(-value(dH1)/(value(R)*value(Tsc[i])) + value(dS1)/value(R))/value(P[i])
    model.K2c[i] == exp(-value(dH2)/(value(R)*value(Tsc[i])) + value(dS2)/value(R))/value(P[i])
    model.K3c[i] == exp(-value(dH3)/(value(R)*value(Tsc[i])) + value(dS3)/value(R))/value(P[i])
    model.K1e[i] == exp(-value(dH1)/(value(R)*value(Tse[i])) + value(dS1)/value(R))/value(P[i])
    model.K2e[i] == exp(-value(dH2)/(value(R)*value(Tse[i])) + value(dS2)/value(R))/value(P[i])
    model.K3e[i] == exp(-value(dH3)/(value(R)*value(Tse[i])) + value(dS3)/value(R))/value(P[i])

# Hydrodynamic Initial Values
for i in model.l:
    model.vg[i] = value(Gas_In_F)*(value(R)*value(Gas_In_T)/value(Gas_In_P))/(0.25*value(pi)*value(Dt)**2)

iv_vmf_a = 1.75/(value(phis)*value(emf)**3)*(value(dp)*value(rhog[1])/value(mug))**2
iv_vmf_b = 150/(value(phis)**2*value(emf)**3)*(1-value(emf))*(value(dp)*value(rhog[1])/value(mug))
iv_vmf_c = -value(dp)**3*value(rhog[1])*(value(rhos)-value(rhog[1]))*value(gc)/value(mug)**2
model.vmf = (-iv_vmf_b + sqrt(iv_vmf_b**2 - 4*iv_vmf_a*iv_vmf_c))/(2*iv_vmf_a)

# =============================================================================
# Solve Model

# -----------------------------------------------------------------------------
# Add arbitary objective function
#model.obj = Objective(expr = 1, sense = maximize)

print('Model Construction Complete     (','%.5f'%(time.time()-splitt),'s)')
splitt = time.time()

# -----------------------------------------------------------------------------
# Set Solver Objects and Options
solver = 'petsc'
opt = SolverFactory(solver)
#opt.options['-snes_monitor']=''
#opt.options['-ksp_monitor']=''
opt.options['-snes_atol']='1e-10'
opt.options['-snes_rtol']='1e-8'

# -----------------------------------------------------------------------------
# Deactivate Constraints

# All q group constraints active

# All r group constraints active through r12
model.eq_r13.deactivate()
model.eq_r14.deactivate()
model.eq_r15.deactivate()
model.eq_r16.deactivate()
model.eq_r17.deactivate()
model.eq_r18.deactivate()
model.eq_r19.deactivate()
model.eq_r20.deactivate()
model.eq_r21.deactivate()
model.eq_r22.deactivate()

model.eq_s1.deactivate()
model.eq_s2.deactivate()
model.eq_s3.deactivate()

model.eq_t1.deactivate()
model.eq_t2.deactivate()

# All group A constraints active

model.eq_b1.deactivate()
model.eq_b2.deactivate()
model.eq_b3.deactivate()
model.eq_b4.deactivate()
model.eq_b5.deactivate()
model.eq_b6.deactivate()
model.eq_b7.deactivate()
model.eq_b8.deactivate()
model.eq_b9.deactivate()
model.eq_b10.deactivate()
model.eq_b11.deactivate()
model.eq_b13.deactivate()
model.eq_b12.deactivate()
model.eq_b14.deactivate()
model.eq_b15.deactivate()
model.eq_b16.deactivate()
model.eq_b17.deactivate()
model.eq_b18.deactivate()

model.eq_c1.deactivate()
model.eq_c2.deactivate()
model.eq_c3.deactivate()
model.eq_c4.deactivate()
model.eq_c5.deactivate()
model.eq_c6.deactivate()
model.eq_c7.deactivate()

model.eq_d1.deactivate()
model.eq_d2.deactivate()
model.eq_d3.deactivate()
model.eq_d4.deactivate()
model.eq_d5.deactivate()
model.eq_d6.deactivate()
model.eq_d7.deactivate()
model.eq_d8.deactivate()
model.eq_d9.deactivate()
model.eq_d10.deactivate()
model.eq_d11.deactivate()

model.eq_e1.deactivate()
model.eq_e2.deactivate()
model.eq_e3.deactivate()
model.eq_e4.deactivate()
model.eq_e5.deactivate()
model.eq_e6.deactivate()

model.eq_f1.deactivate()
model.eq_f2.deactivate()
model.eq_f3.deactivate()
model.eq_f4.deactivate()
model.eq_f5.deactivate()

model.eq_g1.deactivate()

model.eq_h1.deactivate()
model.eq_h2.deactivate()
model.eq_h3.deactivate()
model.eq_h4.deactivate()
model.eq_h5.deactivate()
model.eq_h6.deactivate()
model.eq_h7.deactivate()
model.eq_h8.deactivate()
model.eq_h10.deactivate()
model.eq_h9.deactivate()
model.eq_h11.deactivate()
model.eq_h12.deactivate()

model.eq_i1.deactivate()
model.eq_i2.deactivate()
model.eq_i3.deactivate()
model.eq_i4.deactivate()
model.eq_i5.deactivate()
model.eq_i6.deactivate()

model.eq_j1.deactivate()
model.eq_j2.deactivate()
model.eq_j3.deactivate()
model.eq_j4.deactivate()
model.eq_j5.deactivate()
model.eq_j6.deactivate()
model.eq_j7.deactivate()
model.eq_j8.deactivate()

model.eq_k1.deactivate()
model.eq_k2.deactivate()
model.eq_k3.deactivate()
model.eq_k4.deactivate()
model.eq_k5.deactivate()
model.eq_k6.deactivate()

model.eq_l1.deactivate()

# -----------------------------------------------------------------------------
# Activate IP constraints
model.IP_b5.activate()
model.IP_b6.activate()
model.IP_b9.activate()
model.IP_b10.activate()
model.IP_b15.activate()
model.IP_b16.activate()

# -----------------------------------------------------------------------------
# 1st Initialisation Step - Reactor design and equilibrium solids loading

# Initialise Variables
model.dl = value(Lb)/value(ND)
model.Areact = 0.25*value(pi)*value(Dt)**2
model.lp = value(lhx) + value(dx)
model.Nx = value(Areact)/(value(lp)**2)
model.Ax = value(Areact) - (value(pi)/4)*value(dx)**2*value(Nx)
model.Ao = 1/value(nor)
model.Ahx = value(pi)*value(dx)*value(Lb)*value(Nx)
model.Dte = 4*value(Ax)/(value(pi)*(value(Dt)+value(dx)*value(Nx)))

for i in l:
    model.k1c[i] = value(A1)*value(Tsc[i])*exp(-value(E1)/(value(R)*value(Tsc[i])))
    model.k2c[i] = value(A2)*value(Tsc[i])*exp(-value(E2)/(value(R)*value(Tsc[i])))
    model.k3c[i] = value(A3)*value(Tsc[i])*exp(-value(E3)/(value(R)*value(Tsc[i])))
    model.k1e[i] = value(A1)*value(Tse[i])*exp(-value(E1)/(value(R)*value(Tse[i])))
    model.k2e[i] = value(A2)*value(Tse[i])*exp(-value(E2)/(value(R)*value(Tse[i])))
    model.k3e[i] = value(A3)*value(Tse[i])*exp(-value(E3)/(value(R)*value(Tse[i])))

    model.K1c[i] = exp(-value(dH1)/(value(R)*value(Tsc[i])) + value(dS1)/value(R))/value(P[i])
    model.K2c[i] = exp(-value(dH2)/(value(R)*value(Tsc[i])) + value(dS2)/value(R))/value(P[i])
    model.K3c[i] = exp(-value(dH3)/(value(R)*value(Tsc[i])) + value(dS3)/value(R))/value(P[i])
    model.K1e[i] = exp(-value(dH1)/(value(R)*value(Tse[i])) + value(dS1)/value(R))/value(P[i])
    model.K2e[i] = exp(-value(dH2)/(value(R)*value(Tse[i])) + value(dS2)/value(R))/value(P[i])
    model.K3e[i] = exp(-value(dH3)/(value(R)*value(Tse[i])) + value(dS3)/value(R))/value(P[i])
    
    model.wc['H2O',i] = value(K1c[i])*value(P[i])*value(yc['H2O',i])/value(rhos)
    model.we['H2O',i] = value(K1e[i])*value(P[i])*value(ye['H2O',i])/value(rhos)
    
    # Initial guess of zero Bicarbonate loading to find Carbamate loading
    model.wc['Car',i] = (value(nv)/value(rhos))*sqrt(value(P[i])*value(yc[j,i])*value(K3c[i]))/(1+2*sqrt(value(P[i])*value(yc[j,i])*value(K3c[i])))
    model.we['Car',i] = (value(nv)/value(rhos))*sqrt(value(P[i])*value(ye[j,i])*value(K3e[i]))/(1+2*sqrt(value(P[i])*value(ye[j,i])*value(K3e[i])))

w_a = value(rhos)**2/(value(nv)*value(K2c[1]))
w_b = value(wc['Car',1])*value(rhos)**2/value(K2c[1]) + value(wc['H2O',1])*value(P[1])*value(yc['CO2',1])*value(rhos)**2/value(nv)
w_c = 2*value(wc['H2O',1])*value(wc['Car',i])*value(P[1])*value(yc['CO2',1])*value(rhos)**2/value(nv) - value(wc['H2O',1])*value(rhos)*value(P[1])*value(yc['CO2',1])

for i in l:
    model.wc['Bic',i] = (-w_b + sqrt(w_b**2 - 4*w_a*w_c))/(2*w_a)
    model.we['Bic',i] = (-w_b + sqrt(w_b**2 - 4*w_a*w_c))/(2*w_a)

# Fix Variables
for i in l:
    Gb[i].fixed = True
    P[i].fixed = True
    Phx[i].fixed = True
    Tgb[i].fixed = True
    Tsc[i].fixed = True
    Tse[i].fixed = True
    Thx[i].fixed = True
    
    for j in GasList:
        yb[j,i].fixed = True
        yc[j,i].fixed = True
        ye[j,i].fixed = True

# Solve Step
results = opt.solve(model)

# Report Time and Status
if results.solver.status == SolverStatus.ok:
    print("Initialisation Step  1 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
else:
    print("Initialisation Step  1 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 2nd Initialisation Step - Reaction rate and sorbent flux constraints
if results.solver.status == SolverStatus.ok:
    splitt = time.time()
    
    # Initialise variables
    for i in l:
        model.r1c[i] = 0
        model.r2c[i] = 0
        model.r3c[i] = 0
        model.r1e[i] = 0
        model.r2e[i] = 0
        model.r3e[i] = 0
        for j in GasList:
            model.rgc[j,i] = 0
            model.rge[j,i] = 0
        for j in SolidList:
            model.rsc[j,i] = 0
            model.rse[j,i] = 0

        model.ve[i] = 1.32275*value(vmf)*188*(1.02)*(value(mug)**0.371)/((value(dp)**0.568)*(value(gc)**0.663)*(0.08518*(value(rhos)-value(rhog[i]))+19.09))*(value(Lb)**0.756)/value(Lb)  # Average across bed
        model.dbm[i] = 2.59*(value(gc)**(-0.2))*((value(vg[i])-value(ve[i]))*value(Ax))**0.4
    model.db0 = 1.38*(value(gc)**(-0.2))*((value(vg[1])-value(ve[1]))*value(Ao))**0.4
    model.g1 = 2.56E-2*sqrt(value(Dt)/value(gc))/value(vmf)
    
    db_a = 1+0.3*value(dl)/value(Dt)
    db_b = 0.3*value(dl)*value(g1)/sqrt(value(Dt))
    for i in model.l:
        if i == 1:
            db_c = -0.3*value(dl)*value(dbm[i])/value(Dt) - value(db0)
            
            model.db[i] = ((-db_b + sqrt(db_b**2 - 4*db_a*db_c))/(2*db_a))**2
        else:
            db_c = -0.3*value(dl)*value(dbm[i])/value(Dt) - value(db[i-1])
            
            model.db[i] = ((-db_b + sqrt(db_b**2 - 4*db_a*db_c))/(2*db_a))**2
    
        model.vbr[i] = 0.711*sqrt(value(gc)*value(db[i]))
        model.vb[i] = 1.55*((value(vg[i])-value(vmf))+14.1*(value(db[i])+0.005))*(value(Dte)**0.32)+value(vbr[i])
        model.delta[i] = value(Gas_In_F)*value(V[i])/(value(vb[i])*value(Ax))
        model.Ar[i] = (value(dp)**3)*value(rhog[i])*(value(rhos)-value(rhog[i]))*value(gc)/value(mug)**2
        model.fc[i] = 3*(value(vmf)/value(emf))/(value(vbr[i])-(value(vmf)/value(emf)))
        model.fcw[i] = value(fc[i]) + value(fw)
        model.ed[i] = 1-0.958773*2.05*(1-value(emf))*(value(dp)**0.1)*(value(gc)**0.118)/(2.54*value(mug)**0.066)*value(Lb)**0.043     # Average across bed
        model.e[i] = 1- (1 - value(ed[i]))*(1-value(delta[i]))
        model.Jc[i] = value(fw)*value(delta[i])*value(rhos)*(1-value(ed[i]))*value(vb[i])
        model.Je[i] = value(Jc[i])      # Modify for other boundary conditions
        model.Ksbulk[i] = 0
        for j in SolidList:
            model.Ksbulk_c[j,i] = 0

    # Fix Variables
    for i in l:
        db[i].fixed = True
        delta[i].fixed = True
        ed[i].fixed = True
        fcw[i].fixed = True
        Jc[i].fixed = True

    # Activate Constraints
    model.eq_r13.activate()
    model.eq_r14.activate()
    model.eq_r15.activate()
    model.eq_r16.activate()
    model.eq_r17.activate()
    model.eq_r18.activate()
    model.eq_r19.activate()
    model.eq_r20.activate()
    model.eq_r21.activate()
    model.eq_r22.activate()

    model.eq_b7.activate()
    model.eq_b8.activate()
    model.eq_b13.activate()
    model.eq_b14.activate()
    
    model.eq_f3.activate()

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  2 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  2 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
        
# -----------------------------------------------------------------------------
# 3rd Initialisation Step - Mass balances
if results.solver.status == SolverStatus.ok:
    splitt = time.time()
    
    # Initialise Variables
    for i in l:
        model.cct[i] = 1/value(V[i])
        model.cet[i] = 1/value(V[i])
        model.Kcebs[i] = 3*(1-value(ed[i]))*value(ve[i])/((1-value(delta[i]))*value(ed[i])*value(db[i]))
        for j in GasList:
            model.cb[j,i] = value(Gas_In_y[j])/value(V[i])
            model.cc[j,i] = value(cb[j,i])
            model.ce[j,i] = value(cb[j,i])
            
            model.Kbc[j,i] = 1.32*4.5*(value(vmf)/value(db[i])) + 5.85*((value(D[j,i]))**0.5*(value(gc)**0.25)/(value(db[i])**(5/4)))
            model.Kce[j,i] = 6.77*sqrt(value(ed[i])*value(D[j,i])*value(vbr[i])/value(db[i])**3)
        for j in SolidList:
            model.Kcebs_c[j,i] = value(dl)*value(Ax)*value(delta[i])*value(rhos)*value(Kcebs[i])*(value(wc[j,i])-value(we[j,i]))

    # Fix Variables
    for i in l:
        Kcebs[i].fixed = True
        for j in GasList:
            Kbc[j,i].fixed = True
            Kce[j,i].fixed = True
            Kgbulk_c[j,i].fix(0)
        #for j in SolidList:
        #    Kcebs_c[j,i].fixed = True

    # Activate Constraints
    model.eq_b1.activate()
    model.eq_b2.activate()
    model.eq_b3.activate()
    model.eq_b9.activate()
    model.eq_b10.activate()
    model.eq_b15.activate()
    model.eq_b16.activate()
    
    model.eq_e1.activate()
    model.eq_e2.activate()
    model.eq_e3.activate()
    model.eq_e4.activate()
    model.eq_e5.activate()
    model.eq_e6.activate()
    
    model.eq_f1.activate()
    model.eq_f5.activate()
    
    model.eq_l1.activate()
    
    # Unfix Variables
    for i in l:
        Gb[i].fixed = False
        for j in GasList:
            yb[j,i].fixed = False
            yc[j,i].fixed = False
            ye[j,i].fixed = False
            Kgbulk_c[j,i].fixed = False

    # Deactivate Constraints
    model.IP_b9.deactivate()
    model.IP_b10.deactivate()
    model.IP_b15.deactivate()
    model.IP_b16.deactivate()

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  3 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  3 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 4th Initialisation Step - Enthalpy balances, hydrodynamics and specific enthalpy calculations
if results.solver.status == SolverStatus.ok:
    splitt = time.time()
    
    # Initialise Variables
    model.hs_in = ((value(Solid_In_w['H2O'])+value(Solid_In_w['Bic']))*(value(cpgcse['H2O'])*value(Solid_In_T) + value(dH1)) +
                        value(Solid_In_w['Bic'])*(value(cpgcse['CO2'])*value(Solid_In_T) + value(dH2)) +
                        value(Solid_In_w['Car'])*(value(cpgcse['CO2'])*value(Solid_In_T) + value(dH3))) + value(cps)*value(Solid_In_T)
    model.hhx_in = ((value(HX_In_T)-33.2104)/14170.15 - 0.285)*1e6
    for i in l:
        model.vg[i] = value(Gb[i])*value(V[i])/value(Ax)        
        
        model.hsc[i] = ((value(wc['H2O',i])+value(wc['Bic',i]))*(value(cpgcsc['H2O'])*value(Tsc[i]) + value(dH1)) +
                        value(wc['Bic',i])*(value(cpgcsc['CO2'])*value(Tsc[i]) + value(dH2)) +
                        value(wc['Car',i])*(value(cpgcsc['CO2'])*value(Tsc[i]) + value(dH3))) + value(cps)*value(Tsc[i])
        model.hse[i] = ((value(we['H2O',i])+value(we['Bic',i]))*(value(cpgcse['H2O'])*value(Tse[i]) + value(dH1)) +
                        value(we['Bic',i])*(value(cpgcse['CO2'])*value(Tse[i]) + value(dH2)) +
                        value(we['Car',i])*(value(cpgcse['CO2'])*value(Tse[i]) + value(dH3))) + value(cps)*value(Tse[i])
        model.hhx[i] = ((value(Thx[i])-33.2104)/14170.15 - 0.285)*1e6
        
        model.Hbc[i] = 1.32*4.5*value(vmf)*value(cpg_mol)/(value(db[i])*value(V[i])) +\
                       5.85*sqrt(value(kg)*value(cpg_mol)/value(V[i]))*(value(gc)**0.25)/value(db[i])**(5/4)
        model.Hce[i] = 6.78*sqrt(value(ed[i])*value(vbr[i])*value(kg)*value(cct[i])*value(cpg_mol))/value(db[i])**(3/2)
        
        model.Hgbulk[i] == (6*value(Kd)*value(delta[i])*value(dl)*value(Ax)*(value(cet[i])-(1/value(V[i])))/value(db[i]))*value(cpg_mol)*value(Tgb[i])
        model.Hsbulk[i] = value(Ksbulk[i])*value(hse[i])
        
        model.Ghb[i] = value(Gb[i])*value(Tgb[i])*value(cpg_mol)
        
        model.Jhc[i] = value(Hsbulk[1])/value(Ax)
        model.Jhe[i] = value(Hsbulk[1])/value(Ax)       # Change for different boundary conditions
    
    # Fix Variables
    model.hp.fix(0)
    model.Qhx.fix(0)
    
    for i in l:
        Hbc[i].fixed = True
        Hce[i].fixed = True

    # Activate Constraints
    model.eq_s1.activate()
    model.eq_s2.activate()
    model.eq_s3.activate()
    
    model.eq_t1.activate()
    model.eq_t2.activate()
    
    model.eq_b4.activate()
    model.eq_b5.activate()
    model.eq_b6.activate()
    model.eq_b11.activate()
    model.eq_b12.activate()
    model.eq_b17.activate()
    model.eq_b18.activate()

    model.eq_c1.activate()
    model.eq_c3.activate()

    model.eq_f2.activate()
    model.eq_f4.activate()

    model.eq_h1.activate()
    model.eq_h2.activate()
    model.eq_h3.activate()
    model.eq_h4.activate()
    model.eq_h5.activate()
    model.eq_h6.activate()
    model.eq_h7.activate()
    model.eq_h8.activate()
    model.eq_h10.activate()
    model.eq_h9.activate()
    model.eq_h11.activate()
    model.eq_h12.activate()
    
    # Unfix Variables
    for i in l:
        Tgb[i].fixed = False
        
        model.db[i].fixed = False
        model.ed[i].fixed = False
        model.fcw[i].fixed = False
    
    # Deactivate Constraints
    model.IP_b5.deactivate()
    model.IP_b6.deactivate()

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  4 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  4 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 5th Initialisation Step - Linking hydrodynamics to mass balances, H&M transfer coefficients
if results.solver.status == SolverStatus.ok:
    splitt = time.time()
    
    # Initialise Variables
    for i in l:
        model.delta[i] = value(vg[i])/value(vb[i])
        model.us[i] = value(Je[i])/((1-value(fcw[i])*value(delta[i])-value(delta[i]))*value(rhos)*(1-value(ed[i])))
        
        model.hp[i] = 0.03*value(kg)*((value(ve[i])*value(dp)*value(rhog[i])/value(mug))**1.3)/value(dp)

    # Activate Constraints
    model.eq_c2.activate()
    model.eq_c4.activate()
    model.eq_c6.activate()

    model.eq_i1.activate()
    model.eq_i2.activate()
    model.eq_i3.activate()
    model.eq_i4.activate()
    model.eq_i5.activate()
    model.eq_i6.activate()

    # Unfix Variables
    for i in l:
        model.delta[i].fixed = False
        model.Jc[i].fixed = False
        for j in GasList:
            model.Kbc[j,i].fixed = False
            model.Kce[j,i].fixed = False
        model.Kcebs[i].fixed = False
        model.Hbc[i].fixed = False
        model.Hce[i].fixed = False
        model.hp[i].fixed = False

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  5 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  5 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 6th Initialisation Step - HX Tubes
if results.solver.status == SolverStatus.ok:
    splitt = time.time()
    
    # Initialise Variables
    for i in l:
        model.kpa[i] = (3.58-2.5*value(ed[i]))*value(kg)*((value(kp)/value(kg))**(0.46-0.46*value(ed[i])))
        model.fn[i] = value(vg[i])/value(vmf)
        model.tau[i] = 0.44*((value(dp)*value(gc)/((value(vmf)**2)*((value(fn[i])-value(ah))**2)))**0.14)*((value(dp)/value(dx))**0.225)
        model.fb[i] = 0.33*(((value(vmf)**2)*((value(fn[i])-value(ah))**2)/(value(dp)*value(gc)))**0.14)
        model.hd[i] = 2*sqrt(value(kpa[i])*value(rhos)*value(cps)*(1-value(ed[i])))/sqrt(value(pi)*value(tau[i]))
        model.Pr[i] = value(cpg_mol)*value(mug)/(value(MW[i])*value(kg))
        model.hl[i] = 0.009*(value(Ar[i])**0.5)*(value(Pr[i])**0.33)*value(kg)/value(dp)
        model.ht[i] = value(fb[i])*value(hd[i]) + (1-value(fb[i]))*value(hl[i])

    # Activate Constraints
    model.eq_j1.activate()
    model.eq_j2.activate()
    model.eq_j3.activate()
    model.eq_j4.activate()
    model.eq_j5.activate()
    model.eq_j6.activate()
    model.eq_j7.activate()
    model.eq_j8.activate()
    
    model.eq_k1.activate()
    model.eq_k2.activate()
    model.eq_k3.activate()
    model.eq_k4.activate()
    model.eq_k5.activate()
    model.eq_k6.activate()

    # Unfix Variables
    for i in l:
        Phx[i].fixed = False
        Qhx[i].fixed = False
        Thx[i].fixed = False

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  6 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  6 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 7th Initialisation Step - Linking solids enthalpy and temperature
if results.solver.status == SolverStatus.ok:
    splitt = time.time()

    # Activate Constraints
    model.eq_c5.activate()
    for i in l:
        if i > 1:
            model.eq_c7[i].activate()

    # Unfix Variables
    for i in l:
        Tsc[i].fixed = False
        if i > 1:
            Tse[i].fixed = False

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  7 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  7 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# 8th Initialisation Step - Closing solids enthalpy-temperature loop
if results.solver.status == SolverStatus.ok:
    splitt = time.time()

    # Activate Constraints
    model.eq_c7[1].activate()

    # Unfix Variables
    Tse[1].fixed = False

    # Solve Step
    results = opt.solve(model)

    # Report Time and Status
    if results.solver.status == SolverStatus.ok:
        print("Initialisation Step  8 Complete (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)
    else:
        print("Initialisation Step  8 Failed   (",'%.5f'%(time.time()-splitt),'s)',results.solver.message)

# -----------------------------------------------------------------------------
# Final State Solution - Pressure drop, outlet and diagnostic variables
if results.solver.status == SolverStatus.ok:
    
    # Initialise Variables
    model.Gas_Out_F = value(Gb[value(ND)])
    model.Gas_Out_P = value(P[value(ND)])  
    model.Gas_Out_T = value(Tgb[value(ND)]) 
    for j in GasList:
        model.Gas_Out_y[j] = value(yb[j,value(ND)])  
    
    model.Solid_Out_F = value(Je[value(1)])*value(Ax)       # Change for different boundary conditions
    model.Solid_Out_T = value(Tse[value(1)])
    for j in SolidList:
        model.Solid_Out_w[j] = value(we[j,value(1)])

    model.HX_Out_F = value(HX_In_F)
    model.HX_Out_P = value(Phx[1])
    model.HX_Out_T = value(Thx[1])
    for j in HXList:
        model.HX_Out_y[j] = value(HX_In_y[j])

    # Activate Constraints
    model.eq_d1.activate()
    model.eq_d2.activate()
    model.eq_d3.activate()
    model.eq_d4.activate()
    model.eq_d5.activate()
    model.eq_d6.activate()
    model.eq_d7.activate()
    model.eq_d8.activate()
    model.eq_d9.activate()
    model.eq_d10.activate()
    model.eq_d11.activate()
    
    model.eq_g1.activate()
    
    # Unfix Variables
    for i in l:
        model.P[i].fixed = False
    
    results = opt.solve(model)

# Report Final Time and State
if results.solver.status == SolverStatus.ok:
    print('Solution Completed Successfully (','%.5f'%(time.time()-stime),'s)',results.solver.message)
else:
    print('Solution Failed - Elapsed Time: (','%.5f'%(time.time()-stime),'s)',results.solver.message)
    
print(results)

# -----------------------------------------------------------------------------
# Model Outputs

removal = {}
mbal_tol = {}
mbal_tol['Sorb'] = 1 - Solid_Out_F.value/Solid_In_F.value
for j in GasList:
    if j == 'CO2':
        removal[j] = 1 - Gas_Out_F.value*Gas_Out_y['CO2'].value/(Gas_In_F.value*Gas_In_y['CO2'].value)
        mbal_tol[j] = (Gas_In_F.value*Gas_In_y['CO2'].value + Solid_In_F.value*(Solid_In_w['Car'].value + Solid_In_w['Bic'].value) -  Gas_Out_F.value*Gas_Out_y['CO2'].value - Solid_Out_F.value*(Solid_Out_w['Car'].value + Solid_Out_w['Bic'].value))/(Gas_In_F.value*Gas_In_y['CO2'].value + Solid_In_F.value*(Solid_In_w['Car'].value + Solid_In_w['Bic'].value))
    elif j == 'H2O':
        removal[j] = 1 - Gas_Out_F.value*Gas_Out_y['H2O'].value/(Gas_In_F.value*Gas_In_y['H2O'].value)
        mbal_tol[j] = (Gas_In_F.value*Gas_In_y['H2O'].value + Solid_In_F.value*(Solid_In_w['H2O'].value + Solid_In_w['Bic'].value) -  Gas_Out_F.value*Gas_Out_y['H2O'].value - Solid_Out_F.value*(Solid_Out_w['H2O'].value + Solid_Out_w['Bic'].value))/(Gas_In_F.value*Gas_In_y['H2O'].value + Solid_In_F.value*(Solid_In_w['H2O'].value + Solid_In_w['Bic'].value))
    else:
        mbal_tol[j] = (Gas_In_F.value*Gas_In_y['N2'].value -  Gas_Out_F.value*Gas_Out_y['N2'].value)/(Gas_In_F.value*Gas_In_y['N2'].value)
ebal_gas = value(Gas_In_F)*value(cpg_mol)*value(Gas_In_T) - value(Gas_Out_F)*value(cpg_mol)*value(Gas_Out_T)
ebal_sol = value(Solid_In_F)*value(hs_in) - value(Solid_Out_F)*value(hse[1])        # Modify for boundary conditions
ebal_tol = ebal_gas + ebal_sol + value(Q)

print('Removal:',removal)
print('Mass Balance Tolerance:',mbal_tol)
print('Energy Balance Tolerance:','%.5f'%(ebal_tol),"(Gas:",'%.1f'%(ebal_gas),"Solid:",'%.1f'%(ebal_sol),"HX:",'%.1f'%(value(Q)),")")

Gas_Out_F.display()
Gas_Out_T.display()
Gas_Out_P.display()
Gas_Out_y.display()

Solid_Out_F.display()
Solid_Out_T.display()
Solid_Out_w.display()

HX_Out_F.display()
HX_Out_T.display()
HX_Out_P.display()
HX_Out_y.display()

Q.display()

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:25:32 2016

@author: tburgard

Parameters taken from OxyCombustionMax GAMS code generated during CCSI

"""

from __future__ import division
from __future__ import print_function

__author__ = "Anthony Burgard"
__version__ = "1.0.0"

from pyomo.environ import *

m = ConcreteModel()

# List of Components
comp = ['N2','O2','Ar']

# Critical temperatures
m.Tc = Param(comp,doc='Critical Temperature in K',initialize=\
             {'O2':154.58,'N2':126.2,'Ar':150.86})
m.Pc = Param(comp,doc='Critical Pressure in bar',initialize={'O2':\
    50.45985,'N2':33.943875,'Ar':48.737325})          

Tref = 298.15;      # Reference temperature (K) for Enthalpy calculations
Rsi = 8.314;        # Ideal gas constant (J per K-mol)
m.bigM = 10;        # Big M for 2nd derivative of CEOS
m.EpsilonZ = 1e-6;  # Smallest value for dfdZ in CEOS

# Parameters to help specify the EOS
uEOS = 2
wEOS = -1
omegaA = 0.45724
Inter0 = (uEOS*uEOS - 4*wEOS)**0.5

CpIGTab = {
('1','N2') : 31.12896,
('2','N2') : -1.356E-02,
('3','N2') : 2.678E-05,
('4','N2') : -1.167E-08,
('5','N2') : 0,
('1','O2') : 28.087192,
('2','O2') : -3.678E-06,
('3','O2') : 1.745E-05,
('4','O2') : -1.064E-08,
('5','O2') : 0,
('1','Ar') : 20.790296,
('2','Ar') : -3.209E-05,
('3','Ar') : 5.163E-08,
('4','Ar') : 0,
('5','Ar') : 0,
}
# CpIG units: J/(gmol-K)
m.CpIG = Param(['1','2','3','4','5'],comp,initialize=CpIGTab,doc='Constants\
for Specific heat capacity for ideal gas (from Prop. Gases & Liquids)')  

kEOSTab = {
('N2','N2') : 0,
('O2','N2') : -0.0119,
('Ar','N2') : -0.0026,
('N2','O2') : -0.0119,
('O2','O2') : 0,
('Ar','O2') : 0.0104,
('N2','Ar') : -0.0026,
('O2','Ar') : 0.0104,
('Ar','Ar') : 0,
}
m.kEOS = Param(comp,comp,initialize=kEOSTab,doc='Binary interaction parameters\
 for PR-EOS')

# Pitzer's accentric factor omega (from Prop. Gases & Liquids)
m.omega = Param(comp,doc='Accentric Factor',initialize={'O2':0.021,'N2':0.040,\
'Ar':-0.004})

# Peng-Robinson Equation (PR)
def bEOScalc(m,comp):
    return 0.07780*8.314E-2*m.Tc[comp]/m.Pc[comp]
m.bEOS = Param(comp,initialize=bEOScalc,doc='b for component')

def fwcalc(m,comp):
    return 0.37464 + 1.54226*m.omega[comp] - 0.26992*m.omega[comp]**2
m.fw = Param(comp,initialize=fwcalc,doc='s factor(unique value for each CEOS)')

# Simple function to scale T
def Tscaled(T):
    return T/298.15
    
# Function for ideal gas molar enthalpy in units of kJ/mol
def HIG(x,y,z,T):
    return 0.001 * (x * (Tref**5 * m.CpIG['5','O2'] / 5 * \
                       (Tscaled(T)**5 - 1) + Tref**4 * m.CpIG['4','O2']/ 4 * \
                       (Tscaled(T)**4 - 1) + Tref**3 * m.CpIG['3','O2'] / 3 * \
                       (Tscaled(T)**3 - 1) + Tref**2 * m.CpIG['2','O2'] / 2 * \
                       (Tscaled(T)**2 - 1) + Tref * m.CpIG['1','O2'] * \
                       (Tscaled(T) - 1)) +\
                    y * (Tref**5 * m.CpIG['5','N2'] / 5 * \
                       (Tscaled(T)**5 - 1) + Tref**4 * m.CpIG['4','N2'] / 4 * \
                       (Tscaled(T)**4 - 1) + Tref**3 * m.CpIG['3','N2'] / 3 * \
                       (Tscaled(T)**3 - 1) + Tref**2 * m.CpIG['2','N2'] / 2 * \
                       (Tscaled(T)**2 - 1) + Tref * m.CpIG['1','N2'] * \
                       (Tscaled(T) - 1)) +\
                    z * (Tref**5 * m.CpIG['5','Ar'] / 5 * \
                       (Tscaled(T)**5 - 1) + Tref**4 * m.CpIG['4','Ar'] / 4 * \
                       (Tscaled(T)**4 - 1) + Tref**3 * m.CpIG['3','Ar'] / 3 * \
                       (Tscaled(T)**3 - 1) + Tref**2 * m.CpIG['2','Ar'] / 2 * \
                       (Tscaled(T)**2 - 1) + Tref * m.CpIG['1','Ar'] * \
                       (Tscaled(T) - 1)))

# Function for molar entropy in units of J/(mol K)
def SIG(x,y,z,T):
    return  x*(Tref**3*m.CpIG['4','O2']/3*(Tscaled(T)**3 - 1) +\
    Tref**2 *m.CpIG['3','O2']/2*(Tscaled(T)**2 - 1) + Tref*m.CpIG['2','O2']*\
    (Tscaled(T) - 1) + m.CpIG['1','O2']*log(Tscaled(T))) + \
            y*(Tref**3*m.CpIG['4','N2']/3*(Tscaled(T)**3 - 1) +\
    Tref**2 *m.CpIG['3','N2']/2*(Tscaled(T)**2 - 1) + Tref*m.CpIG['2','N2']*\
    (Tscaled(T) - 1) + m.CpIG['1','N2']*log(Tscaled(T))) + \
            z*(Tref**3*m.CpIG['4','Ar']/3*(Tscaled(T)**3 - 1) +\
    Tref**2 *m.CpIG['3','Ar']/2*(Tscaled(T)**2 - 1) + Tref*m.CpIG['2','Ar']*\
    (Tscaled(T) - 1) + m.CpIG['1','Ar']*log(Tscaled(T)))


def aEOS(j,T):
    return  (1 + m.fw[j]*(1-sqrt(Tscaled(T)*Tref/m.Tc[j])))**2\
                *omegaA*(8.314E-2*m.Tc[j])**2/m.Pc[j]
# Unit of a : bar (m3/kmol)^2

def bmEOS(x,y,z):
    return  x*m.bEOS['O2'] + y*m.bEOS['N2'] + z*m.bEOS['Ar'] 
# Unit of bEOS: m^3/kmol

def amEOS(x,y,z,T):
    return  x*x*sqrt(aEOS('O2',T)*aEOS('O2',T))*(1-m.kEOS['O2','O2']) +\
    x*y*sqrt(aEOS('O2',T)*aEOS('N2',T))*(1-m.kEOS['O2','N2']) +\
    x*z*sqrt(aEOS('O2',T)*aEOS('Ar',T))*(1-m.kEOS['O2','Ar']) +\
    y*x*sqrt(aEOS('N2',T)*aEOS('O2',T))*(1-m.kEOS['N2','O2']) +\
    y*y*sqrt(aEOS('N2',T)*aEOS('N2',T))*(1-m.kEOS['N2','N2']) +\
    y*z*sqrt(aEOS('N2',T)*aEOS('Ar',T))*(1-m.kEOS['N2','Ar']) +\
    z*x*sqrt(aEOS('Ar',T)*aEOS('O2',T))*(1-m.kEOS['Ar','O2']) +\
    z*y*sqrt(aEOS('Ar',T)*aEOS('N2',T))*(1-m.kEOS['Ar','N2']) +\
    z*z*sqrt(aEOS('Ar',T)*aEOS('Ar',T))*(1-m.kEOS['Ar','Ar'])


def bbEOS(x,y,z,T,P):
    return  bmEOS(x,y,z)*(P/101325)/(8.314E-2*Tscaled(T)*Tref) 


def aaEOS(x,y,z,T,P):
    return  amEOS(x,y,z,T)*(P/101325)/(8.314E-2*Tscaled(T)*Tref)**2

def V(x,y,z,T,P,Z):
    return  R*Tscaled(T)*Tref*Z/(P/101325) 


def delta(j,x,y,z,T):
    return  2*sqrt(aEOS(j,T))*( x*sqrt(aEOS('O2',T))*(1 - m.kEOS[j,'O2']) +\
                                y*sqrt(aEOS('N2',T))*(1 - m.kEOS[j,'N2']) +\
                                z*sqrt(aEOS('Ar',T))*(1 - m.kEOS[j,'Ar'])) \
                           / amEOS(x,y,z,T)

def bRatio(j,x,y,z):
    return  m.Tc[j]/m.Pc[j]/(   x*m.Tc['O2']/m.Pc['O2'] +\
                                y*m.Tc['N2']/m.Pc['N2'] +\
                                z*m.Tc['Ar']/m.Pc['Ar']  )

def Inter1(x,y,z,T,P,Z):
    return  Z - bbEOS(x,y,z,T,P)

def Inter2(x,y,z,T,P,Z):
    return  (Z - bbEOS(x,y,z,T,P))/Z


def Inter3(x,y,z,T,P,Z):
    return  (2*Z + bbEOS(x,y,z,T,P)*(uEOS + Inter0)) /\
                 (2*Z+bbEOS(x,y,z,T,P)*(uEOS - Inter0))

def dadT(x,y,z,T):
    return -(Rsi/100/2)*omegaA**0.5*(

    x*x*(1-m.kEOS['O2','O2'])*(m.fw['O2']*(aEOS('O2',T)*m.Tc['O2']/\
    m.Pc['O2'])**0.5 + m.fw['O2']*(aEOS('O2',T)*m.Tc['O2']/m.Pc['O2'])**0.5) +\
    
    x*y*(1-m.kEOS['O2','N2'])*(m.fw['N2']*(aEOS('O2',T)*m.Tc['N2']/\
    m.Pc['N2'])**0.5 + m.fw['O2']*(aEOS('N2',T)*m.Tc['O2']/m.Pc['O2'])**0.5) +\
               
    x*z*(1-m.kEOS['O2','Ar'])*(m.fw['Ar']*(aEOS('O2',T)*m.Tc['Ar']/\
    m.Pc['Ar'])**0.5 + m.fw['O2']*(aEOS('Ar',T)*m.Tc['O2']/m.Pc['O2'])**0.5) +\
               
    y*x*(1-m.kEOS['N2','O2'])*(m.fw['O2']*(aEOS('N2',T)*m.Tc['O2']/\
    m.Pc['O2'])**0.5 + m.fw['N2']*(aEOS('O2',T)*m.Tc['N2']/m.Pc['N2'])**0.5) +\
               
    y*y*(1-m.kEOS['N2','N2'])*(m.fw['N2']*(aEOS('N2',T)*m.Tc['N2']/\
    m.Pc['N2'])**0.5 + m.fw['N2']*(aEOS('N2',T)*m.Tc['N2']/m.Pc['N2'])**0.5) +\
               
    y*z*(1-m.kEOS['N2','Ar'])*(m.fw['Ar']*(aEOS('N2',T)*m.Tc['Ar']/\
    m.Pc['Ar'])**0.5 + m.fw['N2']*(aEOS('Ar',T)*m.Tc['N2']/m.Pc['N2'])**0.5) +\
               
    z*x*(1-m.kEOS['Ar','O2'])*(m.fw['O2']*(aEOS('Ar',T)*m.Tc['O2']/\
    m.Pc['O2'])**0.5 + m.fw['Ar']*(aEOS('O2',T)*m.Tc['Ar']/m.Pc['Ar'])**0.5) +\
               
    z*y*(1-m.kEOS['Ar','N2'])*(m.fw['N2']*(aEOS('Ar',T)*m.Tc['N2']/\
    m.Pc['N2'])**0.5 + m.fw['Ar']*(aEOS('N2',T)*m.Tc['Ar']/m.Pc['Ar'])**0.5) +\
               
    z*z*(1-m.kEOS['Ar','Ar'])*(m.fw['Ar']*(aEOS('Ar',T)*m.Tc['Ar']/\
    m.Pc['Ar'])**0.5 + m.fw['Ar']*(aEOS('Ar',T)*m.Tc['Ar']/m.Pc['Ar'])**0.5))\
         /(Tscaled(T)*Tref)**0.5    
     
def H(x,y,z,T,P,Z):
    return (HIG(x,y,z,T)*bmEOS(x,y,z) +\
        0.1*( Tscaled(T)*Tref*dadT(x,y,z,T)-amEOS(x,y,z,T) )*\
        log(Inter3(x,y,z,T,P,Z))/Inter0 +\
        1E-3*Rsi*Tscaled(T)*Tref*(Z-1)*bmEOS(x,y,z))/bmEOS(x,y,z)
         
def phiEOS(j,x,y,z,T,P,Z):
    return (Inter3(x,y,z,T,P,Z)**(aaEOS(x,y,z,T,P)*(bRatio(j,x,y,z) -\
                delta(j,x,y,z,T)) )*exp(bRatio(j,x,y,z)*( Z - 1 )*\
                bbEOS(x,y,z,T,P)*Inter0))**(1/(bbEOS(x,y,z,T,P)*Inter0))\
                /Inter1(x,y,z,T,P,Z)

def K(j,x,y,z,T,P,Z,x2,y2,z2,T2,P2,Z2):
        return  phiEOS(j,x2,y2,z2,T2,P2,Z2)/\
                phiEOS(j,x,y,z,T,P,Z)

def A_EOS(x,y,z,T,P):
    return -(1 + bbEOS(x,y,z,T,P) - uEOS*bbEOS(x,y,z,T,P))
def B_EOS(x,y,z,T,P):
    return (aaEOS(x,y,z,T,P) + (wEOS-uEOS)*bbEOS(x,y,z,T,P)**2 -\
            uEOS*bbEOS(x,y,z,T,P))
def C_EOS(x,y,z,T,P):
    return -aaEOS(x,y,z,T,P)*bbEOS(x,y,z,T,P) - wEOS*bbEOS(x,y,z,T,P)**2 -\
            wEOS*bbEOS(x,y,z,T,P)**3

def Zequations(m, z1, z2, z3):
    eos = 0 #Use Peng-Robinson EOS
    m.EqZEOSin = Constraint(
        expr =  m.ZEOSin == z1(
        eos, 
        aaEOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P),
        bbEOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P)))
    m.EqZEOSout1 = Constraint(
        expr =  m.ZEOSout1 == z2(
        eos, 
        aaEOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1]),
        bbEOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1])))
    m.EqZEOSout2 = Constraint(
        expr =  m.ZEOSout2 == z3(
        eos, 
        aaEOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2]),
        bbEOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2])))
    
def Zequations2(m):
    m.EqZEOSin = Constraint(expr =  0 == m.ZEOSin**3 +\
    A_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P)*m.ZEOSin**2 +\
    B_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P)*m.ZEOSin +\
    C_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P) )

    m.EqZEOSout1 = Constraint(expr =  0 == m.ZEOSout1**3 +\
    A_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1])*m.ZEOSout1**2 +\
    B_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1])*m.ZEOSout1 +\
    C_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1]) )

    m.EqZEOSout2 = Constraint(expr =  0 == m.ZEOSout2**3 +\
    A_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2])*m.ZEOSout2**2 +\
    B_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2])*m.ZEOSout2 +\
    C_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2]) )

    m.EqZEOSfin = Constraint(expr =  1e-6 <= 3*m.ZEOSin**2 +\
    2*A_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P)*m.ZEOSin +\
    B_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P) )

    m.EqZEOSfout1 = Constraint(expr =  1e-6 <= 3*m.ZEOSout1**2 +\
    2*A_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1])*m.ZEOSout1 +\
    B_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1]) )

    m.EqZEOSfout2 = Constraint(expr =  1e-6 <= 3*m.ZEOSout2**2 +\
    2*A_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2])*m.ZEOSout2 +\
    B_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2]) )
       
    #m.EqVEOSf = Constraint(expr= 6*m.ZEOSin - 2*(1+(1-uEOS)*bbEOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P)) >= 0)        
    m.EqLEOSf = Constraint(expr = 0 >= 6*m.ZEOSin +\
    2*A_EOS(m.In_y['O2'],m.In_y['N2'],m.In_y['Ar'],m.In_T,m.In_P) )
    
    m.EqVEOSflash = Constraint(expr = -10*m.sV <= 6*m.ZEOSout1 +\
    2*A_EOS(m.Out_y[1,'O2'],m.Out_y[1,'N2'],m.Out_y[1,'Ar'],m.Out_T[1],m.Out_P[1]) ) 
    
    m.EqLEOSflash = Constraint(expr = 10*m.sL >= 6*m.ZEOSout2 +\
    2*A_EOS(m.Out_y[2,'O2'],m.Out_y[2,'N2'],m.Out_y[2,'Ar'],m.Out_T[2],m.Out_P[2]) )    

    return
    


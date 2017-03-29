import adolc
from numpy import *
from idaes_models.python_c_wrapper.autogenerate.wrap_decorator import wrapthis

R_gas = 8.314472
TTAG_RHO_VAP = 2

"""
arglist indexes:
0)      T
1)      P
2)      y_CO2
3)      y_H2O
4)      y_N2
5)      y_MEA
6)      y_O2
derivs indexes:
0)      dx/dT
1)      dx/dP
2)      dx/dy_CO2
 
hess indexes:
0)      d2x/dT2
1)      d2x/dTdP
2)      d2x/dP2                                             
3)      d2x/dTdy_CO2
4)      d2x/dPdy_CO2
5)      d2x/dy_CO2dy_CO2
6)      d2x/dy_H2OdT
7)      d2x/dy_H2OdP
8)      d2x/dy_H2Ody_CO2
9)      d2x/dy_H2Ody_H2O
10)...
27) d2x/dy_O2dy_O2          
"""

# The wrapthis decorator indicates to the wrap tool 
# that this method is to be wrapped.
# Internally, the wrap tool checks on the module name of each method,
# if the module name has wrap_decorator in it, it assumes it's a method
# to be wrapped. 
@wrapthis
def cp_ideal(n, at, ra, derivs, hes):
    
    """
    C(i,0), C(i,1)...are constants
    i = chemical component is set of components
    yi = mole fraction of i
    T = temperature
    """
    #TODO: Replace with actual value. 
    TBD = 5.0
    C_i_0 = TBD
    C_i_1 = TBD
    C_i_2 = TBD
    C_i_3 = TBD
    C_i_4 = TBD
    Ci = [C_i_0,C_i_1,C_i_2,C_i_3,C_i_4] # Multiple var names just for clarity, needs to be cleaned up. 
    
    T = ra[at[0]]
    P = ra[at[1]]
    y_CO2 = y0 = ra[at[2]]
    y_H2O = y1 = ra[at[3]]
    y_N2 = y2 = ra[at[4]]
    y_MEA = y3 = ra[at[5]]
    y_O2 = y4 = ra[at[6]] # Multiple var names just for clarity, needs to be cleaned up. 
    yi = [y0,y1,y2,y3,y4]

    cp_ideal = Ci[0]*yi[0] + Ci[1]*yi[1]*T + Ci[2]*yi[2]*T**2 + Ci[3]*yi[3]*T**3

    if derivs:
        derivs[0] = Ci[1]*yi[1] + 2*Ci[2]*yi[2]*T + 3*Ci[3]*yi[3]*T**2
        derivs[1] = Ci[0] + Ci[1]*T + Ci[2]*T**2 + Ci[3]*T**3
        derivs[2:] = [derivs[1]] * (len(derivs) - 2)

        if hes: 
            """
            d**2(cp_ideal)/dT**2 = sum_i((2*C(i,2) + 6*C(i,3)*T)*yi)
            d**2(cp_ideal)/dyi**2 = 0
            d**2(cp_ideal)/dTdyi = C(i,1) + 2*C(i,2)*T + 3*C(i,3)*T**2
            """
            # This could be populated in a smarter way.
            # Column major Hessian matrix indexing + only upper right triangle (including diagonal)
            hes[0] = 2*Ci[2]*yi[2] + 6*Ci[3]*T*yi[3]
            hes[1] = hes[2] = hes [4] = hes[7] = hes[11] = hes[16] = hes[22] = 0 #dPd(anything)
            hes[5] = hes[9] = hes[14] = hes[20] = hes[27] = 0 #../dy**2 [diagonal]
            hes[8] = hes[12] = hes[17] = hes[23] = hes[13] = hes[18] = hes[24] = hes[19] = hes[25] = hes[26] = 0 #../dy**2
            hes[3] = hes[6] = hes[10] = hes[15] = hes[21] = Ci[1] + 2*Ci[2]*T + 3*Ci[3]*T**2 #d**2(cp_ideal)/dTdyi
    return [cp_ideal,derivs,hes]

@wrapthis
def rho_vap (n, at, ra, derivs, hes):
    """
    mi = molecular weight of i (constant in my c file)
    T = temperature
    P = pressure
    R = ideal gas constant (in c file)
    yi = mole fraction of i in vapor
    i is component in set of components
    """
    TBD = 5.0
    m_i_0 = TBD
    m_i_1 = TBD
    m_i_2 = TBD
    m_i_3 = TBD
    m_i_4 = TBD
    mi = [m_i_0,m_i_1,m_i_2,m_i_3,m_i_4]

    T = ra[at[0]]
    P = ra[at[1]]
    y_CO2 = y0 = ra[at[2]]
    y_H2O = y1 = ra[at[3]]
    y_N2 = y2 = ra[at[4]]
    y_MEA = y3 = ra[at[5]]
    y_O2 = y4 = ra[at[6]]

    rho = (m_i_0*y0 + m_i_1*y1 + m_i_2*y2 + m_i_3*y3 + m_i_4*y4)*P/(R_gas*T)
    if derivs:
        derivs[0] = -(m_i_0*y0 + m_i_1*y1 + m_i_2*y2 + m_i_3*y3 + m_i_4*y4)*P/(R*T**2) # d(rho)/dT
        derivs[1] = (m_i_0*y0 + m_i_1*y1 + m_i_2*y2 + m_i_3*y3 + m_i_4*y4)/(R*T) #d(rho)/dP
        for i in range(0,len(derivs)-2):
            derivs[i] = m_list[i]*P/(T*R) # d(rho)/dyi
    if hes:
        # Not sure if this is still column major upper triangle? Maybe I have the matrix wrong, needs to be looked into. 
        hes[0] = 2*(m_i_0*y0 + m_i_1*y1 + m_i_2*y2 + m_i_3*y3 + m_i_4*y4)*P/(R_gas*T**3) # d**2(rho)/dT**2
        hes[2] = 0 # d**2(rho)/dP**2
        hes[5] = hes[9] = hes[14] = hes[20] = hes[27] = 0 #../dy**2 [diagonal] 
        #d**2(rho)/dyidT = -mi*P/(R_gas*T**2)
        #d**2(rho)/dyidP = mi/(R_gas*T)
    return [rho,derivs,hes]
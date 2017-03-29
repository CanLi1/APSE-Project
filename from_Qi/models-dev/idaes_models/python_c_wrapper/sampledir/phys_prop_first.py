from idaes_models.python_c_wrapper.autogenerate.wrap_decorator.wrap_decorator import wrapthis

R_gas = 8.314472

@wrapthis
def V_vap_py_first_five_args(n, at, ra, derivs, hes):
    T = ra[at[0]]
    P = ra[at[1]]
    V = R_gas * T / P
    if derivs:
        derivs[0] = R_gas / P
        derivs[1] = -R_gas * T / (P * P)
        if hes:
            hes[0] = 0
            hes[1] = -R_gas / (P * P)
            hes[2] = 2 * R_gas * T / (P * P * P)
    return [V,derivs,hes]

@wrapthis
def V_vap_py_first_three_args(n, at, ra):
    V = 3
    return [V]

def V_vap_py_first_bad_args(n):
    V = 1
    return [V]
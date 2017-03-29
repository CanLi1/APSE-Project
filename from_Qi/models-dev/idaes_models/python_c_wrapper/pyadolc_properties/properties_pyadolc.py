# This module attempts to implement 2 property methods using py-adolc.
# This is a work in progress and has not been tested but is a proof
# of concept of using pyadolc on Docker.

import adolc
from numpy import *
from idaes_models.python_c_wrapper.autogenerate.wrap_decorator import wrapthis

R_gas = 8.314472
TTAG_RHO_VAP = 2
TTAG_RHO_LIQ = 11

# Component molecular weights in kg/mol
#    [CO2    ,H20    ,N2     ,MEA    ,O2     ]
global mw
mw = [0.04401,0.01802,0.02801,0.06108,0.03200]

#Using pyadolc to migrate rho_vap from C to Python.
@wrapthis
def rho_vap_ad(n, at, ra, derivs, hes):

    # const int n=al->n;
    # static short int tag=TTAG_RHO_VAP;
    # static real r;
    # int i; 
    # adouble T, P, ar;
    # adouble *y = new adouble[n-2]; 
    global mw
    tag = TTAG_RHO_VAP
    r_value = 0 #static real r, return value.
    i = 0
    T = ra[at[0]]
    P = ra[at[1]]
    y_list = np.array([])
    y = adolc.adouble(y_list)
    adolc.trace_on(tag)

    adolc.independent(T) #T <<= al->ra[al->at[0]];
    adolc.independent(P) #P <<= al->ra[al->at[1]];

    # trace_on(tag);
    # T <<= al->ra[al->at[0]];
    # P <<= al->ra[al->at[1]];
    # for(i=2;i<n;++i) y[i-2] <<= al->ra[al->at[i]];
    # ar = 0;
    # for(i=0;i<n-2;++i) ar += mw[i]*y[i];
    # ar *= P/(R_gas*T);
    # ar >>= r;
    # trace_off();
    for i in range(2,n+1):
        y[i-2] = ra[at[i]]
    
    adolc.independent(y)

    # Initial value of r.
    ar = 0
    for i in range(2,n+1):
        ar += mw[i]*y[i]
    ar *= P/(R_gas*T)
    r_value = ar
    # Setting r as dependent variable.
    adolc.dependent(r)
    adolc.trace_off()

    # derivatives call (but without any reuse stuff for now):
    # derivatives(al, tag, 0, NULL, 0);

    # hessian3(tag, n-off, x_cache[tag], &f_cache[tag], g_cache[tag], h_cache[tag]);
    # int hessian3(short tag, int n, double* argument, double *y, double *grad, double** hess) {
    #     //modified version of the adol-c hessian2 driver.
    #     int rc, i,j;
    #     rc = hov_wk_forward(tag,1,n,1,2,n,argument,Xppp[tag],y,Yppp[tag]);
    #     MINDEC(rc,hos_ov_reverse(tag,1,n,1,n,Upp,Zppp[tag]));
    #     for (i=0; i<n; i++){
    #         grad[i] = Yppp[tag][0][i][0];
    #         for (j=0;j<=i;j++) hess[i][j] = Zppp[tag][i][j][1];
    #     }
    #     return rc;
    # }

    # hov_forward(tag,m,n,d,p,x0,X,y0,Y) [C]:
    # hov_wk_forward(1,m,n,d,keep,q,xp,Xppp,yp,Yppp);

    # short int tag; // tape identification
    # int m; // number of dependent variables m
    # int n; // number of independent variables n
    # int d; // highest derivative degree d
    # int p; // number of directions p
    # double x0[n]; // independent vector x0
    # double X[n][p][d]; // tangent matrix X
    # double y0[m]; // dependent vector y0 = F(x0)
    # double Y[m][p][d]; // derivative matrix Y

    # def hov_wk_forward(tape_tag, x, V, keep):
    # (y,W) = hov_wk_forward(tape_tag, x, V, keep)
    # F:R^N -> R^M
    # x is N-vector, y is M-vector
    # D is the order of the derivative
    # V is (N x P x D)-matrix, P directions
    # W is (M x P x D)-matrix, P directional derivatives


    # def hov_ti_reverse(tape_tag, U):


    hes = adolc.hessian_wk_forward(tag,r_value)
    return [r_value,derivs,hes]

@wrapthis
def rho_liq(n, at, ra, derivs, hes):
    global mw
    tag = TTAG_RHO_LIQ
    r_value = 0
    # Skipping reuse. TODO: Add reuse functionality.
    Vst = -1.8218e-6
    a = [0.0, -3.2484e-3, 0.0, -5.35162e-4, 0.0]
    b = [0.0, 1.65, 0.0, -0.451417, 0.0]
    c = [927.0, 793.0, 1.0, 1194.51, 1.0]

    T = ra[at[0]]
    V = rho = 0

    y_list = np.array([])
    adolc.trace_on(tag)

    adolc.independent(T)
    for i in range(1,n): 
        y[i-1] = ra[at[i]];

    for i in range(0,n-1):
        V += y[i]*mw[i]/(a[i]*T*T + b[i]*T + c[i]);

    for i in range(0,n-1):
        rho += y[i]*mw[i];
    rho /= V
    r_value = rho
    adolc.dependent(r_value)
    adolc.trace_off()

    hes = adolc.hessian(tag,r_value)
    return [r_value,derivs,hes]





    
/*
  First pass at a Peng-Robinson equation of state package.  Starting 
  with N2, O2 and Ar for sue with an oxycombustion ASU.

  Author: John Eslick
  Date: November 22, 2016 
  
  Notes:
*/
#include <time.h>
#include <math.h>
#include "funcadd.h"
#undef fprintf
#include <stdio.h>


/***********************************************************************
 * 
 * LOGGING (complie with logging if files defined)
 * 
 **********************************************************************/
#define LOG1 "c.log"

/***********************************************************************
 * 
 * No Real Root Behavior
 * 
 **********************************************************************/
#define FAKE_ROOT
#define FAKE_SMOOTH

/***********************************************************************
 * 
 * CONSTANTS
 * 
 **********************************************************************/

static const double eq_tol=1e-12;
static const double smooth_tol=1e-4;
static const double smooth_eps=1e-4;
static const double one_third=1.0/3.0;
static const double sqrt_3=1.73205080756888;
static const int N_COMP=3; //Number of components

typedef enum{ //Component indexes
    EOS_PR = 0,
    EOS_SRK = 1,
    EOS_END = 2
}eos_indx;

static const double eos_u[EOS_END] = {+2.0, +1.0};
static const double eos_w[EOS_END] = {-1.0, +0.0};

/***********************************************************************
 * 
 * CUBIC ROOT FINDER/DERIVATIVES FUNCTION
 * 
 **********************************************************************/

double curoot(double x){
    //cube root function that can take negatives
    if(x < 0) return -pow(-x, one_third);
    else return pow(x, one_third);
    
}

int cuderiv(eos_indx eos, double A, double B, double z, double *derivs, double *hes){
    const double u = eos_u[eos];
    const double w = eos_w[eos];
    double X, Y, dXdB, dYdA, dYdB;
    
    X = (1 - u)*z*z - (2*w*B - u - 2*u*B)*z + A + 2*w*B + 3*w*B*B;
    Y = 3*z*z - 2*(1 + B - u*B)*z + A + w*B*B - u*B - u*B*B;
    derivs[0] = 0; //derivative with respect to const eos param (not used)
    derivs[1] = (B - z)/Y; //dz/dA
    derivs[2] = X/Y; //dz/dB
    if(hes){
        /* Calculate second derivatives */
        dYdA = 6*z*derivs[1] - 2*(1 + B - u*B)*derivs[1] + 1;
        dYdB = 6*z*derivs[2] - 2*(1 + B - u*B)*derivs[2] - 2*z + 2*u*z + 2*w*B - u - 2*u*B;
        dXdB = 2*(1 - u)*z*derivs[2] - (2*w - 2*u)*z - (2*w*B - u - 2*u*B)*derivs[2] + 2*w + 6*w*B;
        hes[0] = 0; //d2z/deos^2 //eos derivatives should not be used
        hes[1] = 0; //d2z/dAdeos
        hes[3] = 0; //d2z/dBdeos
        hes[2] = -1.0/Y*derivs[1] - (B - z)/(Y*Y)*dYdA; //d2z/dAdA
        hes[4] = 1.0/Y*(1-derivs[2]) - (B - z)/(Y*Y)*dYdB; //d2z/dAdB
        hes[5] = 1.0/Y*dXdB - X/(Y*Y)*dYdB; //d2z/dBdB
    }
    return 0;
}

double cubic_root(int phase, eos_indx eos, double A, double B, 
                  double *derivs, double *hes){
    /* 
     * Find solution to the cubic equation:
     * 0 = z^3 - (1+B-uB)z^2 + (A+wB-uB-uB^2)z - AB-wB^2 - wB^3
     * 
     * Also using the definition below (no a, a = 1)
     * 
     * 0 = z^3 + b*z2 + c*z + d
     * 
     * Return either what should be the liquid root if phase == 0 or the
     * vapor root if phase == 1, and the derivatives of the returned 
     * root
     * 
     * The root finding approach given in CRC Standard Mathematical 
     *  Tables and Formulae 32nd edition pg 68 was used to identify the 
     *  roots.
     * 
     */
    int i;
    const double u = eos_u[eos];
    const double w = eos_w[eos];
    static double b, c, d;
    static double F, G, H, I, J, K, M, N, P, R, S, T, U;
    static double dPdb, dPdc, dPdd, dFdb, dFdc, dFdd, dGdb, dGdc, dGdd;
    static double dHdb, dHdc, dHdd, dRdb, dRdc, dRdd, dTdb, dTdc, dTdd;
    static double dSdb, dSdc, dSdd, dUdb, dUdc, dUdd;
    static double dzdb, dzdc, dzdd, dbdA, dbdB, dcdA, dcdB, dddA, dddB;
    static double cim;
    static double d2bAA, d2cAA, d2dAA, d2bAB, d2cAB, d2dAB;
    static double d2bBB, d2cBB, d2dBB;
    static double d2Hbb, d2Hbc, d2Hbd, d2Hcc, d2Hcd, d2Hdd;
    static double d2Gbb, d2Gbc, d2Gbd, d2Gcc, d2Gcd, d2Gdd;
    static double d2Fbb, d2Fbc, d2Fbd, d2Fcc, d2Fcd, d2Fdd;
    static double d2Rbb, d2Rbc, d2Rbd, d2Rcc, d2Rcd, d2Rdd;
    static double d2Tbb, d2Tbc, d2Tbd, d2Tcc, d2Tcd, d2Tdd;
    static double d2Sbb, d2Sbc, d2Sbd, d2Scc, d2Scd, d2Sdd;
    static double d2Ubb, d2Ubc, d2Ubd, d2Ucc, d2Ucd, d2Udd;
    static double d2zbb, d2zbc, d2zbd, d2zcc, d2zcd, d2zdd;
    
    double z, z2 = 0, zr[3], zi[3]={0.0, 0.0, 0.0};
    char fake=0;
    double inflect;
    double resid;
    double smooth_loc;
    #ifdef LOG1
        FILE *file_ptr;
    #endif
    
    b = -(1 + B - u*B);
    c = A + w*B*B - u*B - u*B*B;
    d = -A*B - w*B*B - w*B*B*B;
    inflect = -b/3.0;
    F = (3*c - b*b)/3.0;
    G = (2*b*b*b - 9*b*c + 27*d)/27.0;
    H = G*G/4.0 + F*F*F/27.0;
    P = -b/3.0;
    
    if(H <= 0.0){
        //Three roots, can include double or triple roots too
        I = sqrt(G*G/4.0 - H);
        J = curoot(I);
        K = acos(-G/2.0/I);
        M = cos(K/3);
        N = sqrt_3*sin(K/3.0);
        zr[0] = P + 2*J*M;
        zr[1] = P - J*(M + N);
        zr[2] = P - J*(M - N);
        //Sort lowest to highest, may be a way to not need to sort and
        //know which order they should be in espessially give a = 1, but
        //I don't want to think about that now.
        if(zr[0] > zr[1]){z = zr[1]; zr[1] = zr[0]; zr[0] = z;}
        if(zr[1] > zr[2]){z = zr[2]; zr[2] = zr[1]; zr[1] = z;}
        if(zr[0] > zr[1]){z = zr[1]; zr[1] = zr[0]; zr[0] = z;}
        //Vapor is the highest and liquid is the lowest
        if(phase == 0) z = zr[0]; //liquid
        else if(phase == 1) z = zr[2]; //vapor
        else z = 0.0/0.0;
    }
    else{
        R = -G/2.0 + sqrt(H);
        S = curoot(R);
        T = -G/2.0 - sqrt(H);
        U = -curoot(-T);
        zr[0] = S + U + P;
        zr[1] = -(S + U)/2.0 + P;
        zr[2] = -(S + U)/2.0 + P;
        zi[1] = -sqrt_3*(S - U)/2.0;
        zi[2] = +sqrt_3*(S - U)/2.0;
        if(fabs(zr[0] - inflect) < eq_tol) z = zr[0];//liquid/vapor almost same
        else if(phase == 0 && zr[0] < inflect) z = zr[0];//liquid root
        else if(phase == 1 && zr[0] > inflect) z = zr[0];//vapor root
        else{
            #ifndef FAKE_ROOT
                z = 0.0/0.0;
            #else
                if(phase == 0){ z = zr[2]; z2 = zi[2];}
                else { z = zr[1]; z2 = zi[1];}
                fake = 1;
            #endif
        }
    }
    if(fake) z = z + z2;
    resid = z*z*z + b*z*z + c*z + d; // Claculate the z residul
    /* Derivatives */
    /* Calculate first derivatives with respect to coefficients A and B
     * also fill 0 for derivatives with respect to the constant EOS
     * parameter */
    if(derivs && z == z && !fake) cuderiv(eos, A, B, z, derivs, hes);
    else if(derivs && fake){
        //Didn't move this to a seperate function because it uses
        //a ton of stuff already calculated above.
        dPdb = -one_third;
        dPdc = 0.0;
        dPdd = 0.0;
        dFdb = -2.0*b/3;
        dFdc = 1.0;
        dFdd = 0.0;
        dGdb = (6.0*b*b - 9*c)/27.0;
        dGdc = -9*b/27.0;
        dGdd = 1.0;
        dHdb = 0.5*G*dGdb + F*F*dFdb/9.0;
        dHdc = 0.5*G*dGdc + F*F*dFdc/9.0;
        dHdd = 0.5*G*dGdd + F*F*dFdd/9.0;
        dRdb = -0.5*dGdb + 0.5/sqrt(H)*dHdb;
        dRdc = -0.5*dGdc + 0.5/sqrt(H)*dHdc;
        dRdd = -0.5*dGdd + 0.5/sqrt(H)*dHdd;
        dTdb = -0.5*dGdb - 0.5/sqrt(H)*dHdb;
        dTdc = -0.5*dGdc - 0.5/sqrt(H)*dHdc;
        dTdd = -0.5*dGdd - 0.5/sqrt(H)*dHdd;
        dSdb = one_third/pow(R*R, one_third)*dRdb;
        dSdc = one_third/pow(R*R, one_third)*dRdc;
        dSdd = one_third/pow(R*R, one_third)*dRdd;
        dUdb = one_third/pow(T*T, one_third)*dTdb;
        dUdc = one_third/pow(T*T, one_third)*dTdc;
        dUdd = one_third/pow(T*T, one_third)*dTdd;
        if(phase==0) cim = 1;
        else cim = -1;
        dzdb = -0.5*dSdb - 0.5*dUdb + dPdb + cim*sqrt_3/2.0*(dSdb - dUdb);
        dzdc = -0.5*dSdc - 0.5*dUdc + dPdc + cim*sqrt_3/2.0*(dSdc - dUdc);
        dzdd = -0.5*dSdd - 0.5*dUdd + dPdd + cim*sqrt_3/2.0*(dSdd - dUdd);
        dbdA = 0.0;
        dbdB = u - 1;
        dcdA = 1.0;
        dcdB = 2.0*w*B - u -2*u*B;
        dddA = -B;
        dddB = -A - 2*w*B - 3*w*B*B;
        derivs[0] = 0.0;
        derivs[1] = dzdb*dbdA + dzdc*dcdA + dzdd*dddA;
        derivs[2] = dzdb*dbdB + dzdc*dcdB + dzdd*dddB;
        if(hes){
            d2bAA = 0.0;
            d2cAA = 0.0;
            d2dAA = 0.0;
            d2bAB = 0.0;
            d2bAB = 0.0;
            d2dAB = -1.0;
            d2bBB = 0.0;
            d2cBB = 2*w - 2*u;
            d2dBB = -2*w - 6*w*B;
            d2Fbb = -2.0/3.0;
            d2Fbc = 0;
            d2Fbd = 0;
            d2Fcc = 0;
            d2Fcd = 0;
            d2Fdd = 0;
            d2Gbb = 12*b/27.0;
            d2Gbc = -one_third;
            d2Gbd = 0;
            d2Gcc = 0;
            d2Gcd = 0;
            d2Gdd = 0;
            d2Hbb = 0.5*dGdb*dGdb + 0.5*G*d2Gbb + 2.0/9.0*F*dFdb*dFdb + F*F/9.0*d2Fbb;
            d2Hbc = 0.5*dGdb*dGdc + 0.5*G*d2Gbc + 2.0/9.0*F*dFdb*dFdc + F*F/9.0*d2Fbc;
            d2Hbd = 0.5*dGdb*dGdd + 0.5*G*d2Gbd + 2.0/9.0*F*dFdb*dFdd + F*F/9.0*d2Fbd;
            d2Hcc = 0.5*dGdc*dGdc + 0.5*G*d2Gcc + 2.0/9.0*F*dFdc*dFdc + F*F/9.0*d2Fcc;
            d2Hcd = 0.5*dGdc*dGdd + 0.5*G*d2Gcd + 2.0/9.0*F*dFdc*dFdd + F*F/9.0*d2Fcd;
            d2Hdd = 0.5*dGdd*dGdd + 0.5*G*d2Gdd + 2.0/9.0*F*dFdd*dFdd + F*F/9.0*d2Fdd;
            d2Rbb = -0.5*d2Gbb + 0.5/sqrt(H)*d2Hbb - 0.25*pow(H, -2.0/3.0)*dHdb*dHdb;
            d2Rbc = -0.5*d2Gbc + 0.5/sqrt(H)*d2Hbc - 0.25*pow(H, -2.0/3.0)*dHdb*dHdc;
            d2Rbd = -0.5*d2Gbd + 0.5/sqrt(H)*d2Hbd - 0.25*pow(H, -2.0/3.0)*dHdb*dHdd;
            d2Rcc = -0.5*d2Gcc + 0.5/sqrt(H)*d2Hcc - 0.25*pow(H, -2.0/3.0)*dHdc*dHdc;
            d2Rcd = -0.5*d2Gcd + 0.5/sqrt(H)*d2Hcd - 0.25*pow(H, -2.0/3.0)*dHdc*dHdd;
            d2Rdd = -0.5*d2Gdd + 0.5/sqrt(H)*d2Hdd - 0.25*pow(H, -2.0/3.0)*dHdd*dHdd;
            d2Tbb = -0.5*d2Gbb - 0.5/sqrt(H)*d2Hbb + 0.25*pow(H, -2.0/3.0)*dHdb*dHdb;
            d2Tbc = -0.5*d2Gbc - 0.5/sqrt(H)*d2Hbc + 0.25*pow(H, -2.0/3.0)*dHdb*dHdc;
            d2Tbd = -0.5*d2Gbd - 0.5/sqrt(H)*d2Hbd + 0.25*pow(H, -2.0/3.0)*dHdb*dHdd;
            d2Tcc = -0.5*d2Gcc - 0.5/sqrt(H)*d2Hcc + 0.25*pow(H, -2.0/3.0)*dHdc*dHdc;
            d2Tcd = -0.5*d2Gcd - 0.5/sqrt(H)*d2Hcd + 0.25*pow(H, -2.0/3.0)*dHdc*dHdd;
            d2Tdd = -0.5*d2Gdd - 0.5/sqrt(H)*d2Hdd + 0.25*pow(H, -2.0/3.0)*dHdd*dHdd;
            d2Rbb = -2.0/9.0/curoot(R*R*R*R*R)*dRdb*dRdb + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rbb;
            d2Rbc = -2.0/9.0/curoot(R*R*R*R*R)*dRdb*dRdc + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rbc;
            d2Rbd = -2.0/9.0/curoot(R*R*R*R*R)*dRdb*dRdd + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rbd;
            d2Rcc = -2.0/9.0/curoot(R*R*R*R*R)*dRdc*dRdc + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rcc;
            d2Rcd = -2.0/9.0/curoot(R*R*R*R*R)*dRdc*dRdd + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rcd;
            d2Rdd = -2.0/9.0/curoot(R*R*R*R*R)*dRdd*dRdd + 1.0/3.0*pow(R*R, -2.0/3.0)*d2Rdd;
            d2Tbb = -2.0/9.0/curoot(-T*T*T*T*T)*dTdb*dRdb + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tbb;
            d2Tbc = -2.0/9.0/curoot(-T*T*T*T*T)*dTdb*dRdc + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tbc;
            d2Tbd = -2.0/9.0/curoot(-T*T*T*T*T)*dTdb*dRdd + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tbd;
            d2Tcc = -2.0/9.0/curoot(-T*T*T*T*T)*dTdc*dRdc + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tcc;
            d2Tcd = -2.0/9.0/curoot(-T*T*T*T*T)*dTdc*dRdd + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tcd;
            d2Tdd = -2.0/9.0/curoot(-T*T*T*T*T)*dTdd*dRdd + 1.0/3.0*pow(T*T, -2.0/3.0)*d2Tdd;
            d2zbb = -0.5*d2Sbb - 0.5*d2Ubb + cim*sqrt_3/2.0*(d2Sbb + d2Ubb);
            d2zbc = -0.5*d2Sbc - 0.5*d2Ubc + cim*sqrt_3/2.0*(d2Sbc + d2Ubc);
            d2zbd = -0.5*d2Sbd - 0.5*d2Ubd + cim*sqrt_3/2.0*(d2Sbd + d2Ubd);
            d2zcc = -0.5*d2Scc - 0.5*d2Ucc + cim*sqrt_3/2.0*(d2Scc + d2Ucc);
            d2zcd = -0.5*d2Scd - 0.5*d2Ucd + cim*sqrt_3/2.0*(d2Scd + d2Ucd);
            d2zdd = -0.5*d2Sdd - 0.5*d2Udd + cim*sqrt_3/2.0*(d2Sdd + d2Udd);
        }
    }
    
    #ifdef LOG1
    //Some stuff to check input/output when used in solver.
    file_ptr = fopen(LOG1, "a");
    fprintf(file_ptr, 
        "phase, A, B, z, z0, z1, zc2, resid = %d, %f, %f, %f, (%f, %f), (%f, %f), (%f, %f), %f\n", 
        phase, A, B, z, zr[0], zi[0], zr[1], zi[1], zr[2], zi[2], resid);
    fclose(file_ptr);
    #endif
    return z;
}

/***********************************************************************
 * 
 * LIQUID/VAPOR ROOT PROPERTY FUNCTIONS
 * 
 **********************************************************************/

real ceos_z_liq(arglist *al){
    /* This finds the liquid root for a cubic EOS given the coefficents of
     * 0 = z^3 + b*z^2 + c*z + d
     * 
     * If also calculates first and second derivative.  Returns 0.0/0.0
     * in the case the liquid root does not exist.  Can be used as a
     * building block for cubic EOS property methods
     */
    double A, B;
    eos_indx eos;
    
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    
    return cubic_root(0, eos, A, B, al->derivs, al->hes);
}

real ceos_z_vap(arglist *al){
    /* This finds the vapor root for a cubic EOS given the coefficents of
     * 0 = z^3 + b*z^2 + c*z + d
     * 
     * If also calculates first and second derivative.  Returns 0.0/0.0
     * in the case the liquid root does not exist.  Can be used as a
     * building block for cubic EOS property methods
     */
    double A, B;
    eos_indx eos;
    
    eos = (eos_indx)(al->ra[al->at[0]] + 0.5);  //cast to int is floor
    A = al->ra[al->at[1]];
    B = al->ra[al->at[2]];
    
    return cubic_root(1, eos, A, B, al->derivs, al->hes);
}

void funcadd(AmplExports *ae){
    /* Arguments for addfunc (this is not fully detailed see funcadd.h)
     * 1) Name of function in AMPL
     * 2) Function pointer to C function
     * 3) see FUNCADD_TYPE enum in funcadd.h
     * 4) Number of arguments (the -1 is variable arg list length)
     * 5) Void pointer to function info
     */
    //Real value function, and suppress anoying warings that, I think 
    //happen when this gets called on an already loaded library.
    int t = FUNCADD_REAL_VALUED;  //(FUNCADD_NO_DUPWARN doesn't work)
    addfunc("ceos_z_vap", (rfunc)ceos_z_vap, t, -1, NULL);
    addfunc("ceos_z_liq", (rfunc)ceos_z_liq, t, -1, NULL);
}











//#define function function_old
//#include "adolc/adolc.h"
//#undef function
#define No_AE_redefs
#include "funcadd.h"
#define Printf printf

#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void *dlhandle; // dynamic library handle

/* Function pointer type for AMPL user function */ 
typedef real (*auf_t)(arglist *al);
typedef real (*croot_t)(int phase, int eos, double A, double B, 
                        double *derivs, double *hes);


double residual(int eos, double A, double B, double z){
    double b,c,d,u,w;
    if(eos==0){
        u = 2;
        w = -1;}
    else{
        u = 1;
        w = 0;}
    b = -(1 + B - u*B);
    c = A + w*B*B - u*B - u*B*B;
    d = -A*B - w*B*B - w*B*B*B;
    return z*z*z + b*z*z + c*z + d;
}
                        
int main(int argc, char **argv){
    char *filename; // .so file name
    croot_t f;
    double z, z2, z3;
    double derivs[3], derivs2[3], derivs3[3];
    double fd_dzdA, fd_dzdB;
    double A, B, AA, BB, T, P;
    int eos, i;
    double h=1e-6;
    double hes[6], hes2[6];
    int phase;

    phase = 0;
    eos = 0;
    A = 0.10;
    B = 0.014;
    
    /* Print some heading */
    printf("\n\n");
    printf("+-----------------------------+\n");
    printf("|   Property library tests    |\n");
    printf("+-----------------------------+\n\n");
    /* Check that file may have been specified */
    if(argc < 2){
        printf("Must specify library file name.\n");
        return 1;
    }
    /* Load so file */
    filename = argv[1];
    printf("Library file: %s\n", filename);
    dlhandle = dlopen(filename, RTLD_NOW|RTLD_GLOBAL);
    if(!dlhandle){ // if didn't load print error
        printf("%s\n", dlerror());
        return 2;
    }
    printf("Loaded library\n");
    
    f = (croot_t)dlsym(dlhandle, "cubic_root");
    
    z = f(phase, eos, A, B, derivs, hes);
    printf("phase: %d\n", phase);
    printf("eos: %d, 0 = PR, 1 = SRK\n", eos);
    printf("A, B: %f, %f\n", A, B);
    printf("root: %f\n", z);
    printf("residual: %e\n", residual(eos, A, B, z));
    printf("\nFirst Derivatives\n");
    printf("dz/dA: %f\n", derivs[1]); 
    printf("dz/dB: %f\n", derivs[2]);
    printf("\nSecond Derivatives\n");
    printf("d2z/dAdA: %f\n", hes[2]);
    printf("d2z/dAdB: %f\n", hes[4]);
    printf("d2z/dBdB: %f\n", hes[5]);
    
    /* Finite difference check */
    printf("\nCentered finite difference approx. for first derivatives\n");
    z2 = f(phase, eos, A+h/2.0, B, derivs, hes2);
    z3 = f(phase, eos, A-h/2.0, B, derivs, hes2);
    fd_dzdA = (z2 - z3)/h;
    printf("dz/dA: %f\n", fd_dzdA); 
    
    z2 = f(phase, eos, A, B+h/2.0, derivs, hes2);
    z3 = f(phase, eos, A, B-h/2.0, derivs, hes2);
    fd_dzdB = (z2 - z3)/h;
    printf("dz/dB: %f\n", fd_dzdB);
        
    /* Finite difference check hessian*/
    printf("\nHessian Check\n");
    printf("(Exact, Centered finite difference of exact first derivatives)\n");
    z2 = f(phase, eos, A+h/2.0, B, derivs2, hes2);
    z3 = f(phase, eos, A-h/2.0, B, derivs3, hes2);

    printf("d2z/dAdA = %f\n", (derivs2[1] - derivs3[1])/h);
    printf("d2z/dAdB = %f\n", (derivs2[2] - derivs3[2])/h);

    z2 = f(phase, eos, A, B+h/2.0, derivs2, hes2);
    z3 = f(phase, eos, A, B-h/2.0, derivs3, hes2);
    
    printf("d2z/dBdB = %f\n", (derivs2[2] - derivs3[2])/h);
    printf("\n");
    
    for(i=0; i<200; ++i){
        AA = A - 0.0001*i;
        BB = B;
        z = f(phase, eos, AA, BB, derivs, hes);
        printf("%d, %f, %f, %f, %f, %f, %f\n", i, AA, BB, z,  derivs[1], derivs[2], residual(eos, AA, BB, z));
    }

    dlclose(dlhandle);
    return 0;
}




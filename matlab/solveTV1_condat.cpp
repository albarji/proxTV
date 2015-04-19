#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveTV1_condat.cpp

   Solves the TV-L1 proximity problem by applying Condat's segment construction method.

   Parameters:
     - 0: reference signal y.
     - 1: lambda penalty.
     
   Outputs:
     - 0: primal solution x.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x=NULL,*y;
    float lambda;
    int M,N,nn,i;
    
    #define FREE \
        if(!nlhs) free(x);
    
    #define CANCEL(txt) \
        printf("Error in solveTV1_condat: %s\n",txt); \
        if(x) free(x); \
        return;
    
    /* Check input correctness */
    if(nrhs < 2){CANCEL("not enought inputs");}
    if(!mxIsClass(prhs[0],"double")) {CANCEL("input signal must be in double format")}

    /* Create output arrays */
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    nn = (M > N) ? M : N;
    if(nlhs >= 1){
        plhs[0] = mxCreateDoubleMatrix(nn,1,mxREAL);
        x = mxGetPr(plhs[0]);
    }
    else x = (double*)malloc(sizeof(double)*nn);

    /* Retrieve input data */
    y = mxGetPr(prhs[0]);
    lambda = mxGetScalar(prhs[1]);
    
    /* Run Condat's method  */
    TV1D_denoise(y, x, nn, lambda);
    
    /* Free resources */
    FREE
    
    return;
    
    #undef FREE
    #undef CANCEL
}



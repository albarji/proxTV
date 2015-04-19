#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveTVp_GPFW.cpp

   Solves the general TV-Lp proximity problem by applying a Gradient Projection + Frank-Wolfe method.

   Parameters:
     - 0: reference signal y.
     - 1: lambda penalty.
     - 2: p norm.
     
   Outputs:
     - 0: primal solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: dual gap.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x=NULL,*y,*info=NULL;
    double lambda,p;
    int M,N,nn,i;
    
    #define FREE \
        if(!nlhs) free(x);
    
    #define CANCEL(txt) \
        printf("Error in solveTVp_GP: %s\n",txt); \
        if(x) free(x); \
        if(info) free(info); \
        return;
    
    /* Check input correctness */
    if(nrhs < 3){CANCEL("not enought inputs");}
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
    if(nlhs >= 2){
        plhs[1] = mxCreateDoubleMatrix(N_INFO,1,mxREAL);
        info = mxGetPr(plhs[1]);
    }

    /* Retrieve input data */
    y = mxGetPr(prhs[0]);
    lambda = mxGetScalar(prhs[1]);
    p = mxGetScalar(prhs[2]);
    
    /* Run GP+FW method  */
    GPFW_TVp(y, lambda, x, info, nn, p, NULL);
    
    /* Free resources */
    FREE
    
    return;
}


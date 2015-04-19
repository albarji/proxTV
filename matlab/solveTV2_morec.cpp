#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveTV2_morec.cpp

   Solves the TV-L2 proximity problem by applying a More-Sorensen algorithm.

   Parameters:
     - 0: reference signal y.
     - 1: lambda penalty.
     
   Outputs:
     - 0: primal solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: dual gap.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x,*y,*info;
    double lambda;
    int M,N,nn,i;

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
    
    /* Run Projected Newton */
    more_TV2(y,lambda,x,info,nn);
    
    /* Free resources */
    if(!nlhs) free(x);
    
    return;
}



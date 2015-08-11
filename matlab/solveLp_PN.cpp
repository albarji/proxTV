#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveLp_PN.cpp

   Matlab interface to Lp proximity solver using Projected Newton.

   Parameters:
     - 0: reference signal y.
     - 1: lambda penalty.
     - 2: p degree of the norm.
     - 3: (optional) maximum number of iterations in the solver.
        If not specified the algorithm will run until a certain
        dual gap value is attained
     
   Outputs:
     - 0: primal solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: stopping tolerance.
     - 2: array with dual gap at each iteration. Only returned if
        a maximum number of iterations was provided.
     - 3: array with computation times in microseconds at each iteration.
        Only returned if a maximum number of iterations was provided.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x=NULL,*y,*info=NULL,*gaps=NULL,*times=NULL;
    float lambda,p;
    int nds,M,N,nn,i,maxIters;
    
    #define FREE \
        if(!nlhs) free(x);
    
    #define CANCEL(txt) \
        printf("Error in solveLp_PN: %s\n",txt); \
        if(x) free(x); \
        return;
        
    /* Check input correctness */
    if(nrhs < 3){CANCEL("not enought inputs");}
    if(!mxIsClass(prhs[0],"double")) {CANCEL("input signal must be in double format")}
    /* Get number of dimensions of input signal */
    nds = mxGetNumberOfDimensions(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    if(nds > 2 || (M > 1 && N > 1)){CANCEL("Lp proximity only supported for vector data")}

    /* Retrieve input data */
    y = mxGetPr(prhs[0]);
    lambda = mxGetScalar(prhs[1]);
    p = mxGetScalar(prhs[2]);
    if(nrhs >= 4)
        maxIters = (int)mxGetScalar(prhs[3]);
    else
        maxIters = 0;
        
    /* Create output arrays */
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
    else info = NULL;
    if(nlhs >= 3 && maxIters > 0){
        plhs[2] = mxCreateDoubleMatrix(maxIters+1,1,mxREAL);
        gaps = mxGetPr(plhs[2]);
    }
    if(nlhs >= 4 && maxIters > 0){
        plhs[3] = mxCreateDoubleMatrix(maxIters+1,1,mxREAL);
        times = mxGetPr(plhs[3]);
    }
    
    /* Run solver */
    if ( maxIters > 0 )
        PN_LPp(y, lambda, x, info, nn, p, NULL, 0, 0, maxIters, gaps, times);
    else
        PN_LPp(y, lambda, x, info, nn, p, NULL, 0);
    
    /* Free resources */
    FREE
    
    return;
}



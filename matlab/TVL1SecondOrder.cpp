#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* TVL1SecondOrder.cpp

   Matlab interface to 2nd order TV-L1 proximity solver using Projected Newton.

   Parameters:
     - 0: reference signal y.
     - 1: lambda weight vector (length n-2)
     
   Outputs:
     - 0: primal solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: stopping tolerance.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x=NULL,*y,*info=NULL;
    double* lambda; double p;
    int nds,M,N,nn,i;
    
    #define FREE \
        if(!nlhs) free(x);
    
    #define CANCEL(txt) \
        printf("Error in TVL1Weighted: %s\n",txt); \
        if(x) free(x); \
        return;
        
    /* Check input correctness */
    if(nrhs < 2){CANCEL("not enought inputs");}
    if(!mxIsClass(prhs[0],"double")) {CANCEL("input signal must be in double format")}
    /* Get number of dimensions of input signal */
    nds = mxGetNumberOfDimensions(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    if(nds > 2 || (M > 1 && N > 1)){CANCEL("Lp proximity only supported for vector data")}

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

    /* Retrieve input data */
    y = mxGetPr(prhs[0]);
    lambda = mxGetPr(prhs[1]);
    
    /* Run solver */
    PN_TV1_Trend2_Weighted(y, lambda, x, info, nn, SIGMA, NULL);
    
    /* Free resources */
    FREE
    
    return;
}



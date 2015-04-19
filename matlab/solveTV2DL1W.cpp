#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveTV2DL1W_tautString.cpp

   Matlab interface to weighted 2D-TV-L1 proximity solver using taut string method.

   Parameters:
     - 0: reference signal y of size MxN
     - 1: regularization weight along cols, as a matrix of size (M-1)xN
     - 2: regularization weight along rows, as a matrix of size Mx(N-1)
     - 3: (optional) number of cores to use (default: as defined by environment variable OMP_NUM_THREADS)
     - 4: (optional) maximum number of iterations to run (default: as defined in TVopt.h)
     
   Outputs:
     - 0: primal solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: stopping tolerance.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    double *x=NULL,*y, *W1, *W2, *info=NULL;
    int nds,M,N,nn,i,ncores,maxIters;
    
    #define FREE \
        if(!nlhs) free(x);
    
    #define CANCEL(txt) \
        printf("Error in TVLsolveTV2DL1W_tautString: %s\n",txt); \
        if(x) free(x); \
        return;
        
    /* Check input correctness */
    if(nrhs < 3){CANCEL("not enought inputs");}
    if(!mxIsClass(prhs[0],"double")) {CANCEL("input signal must be in double format")}
    /* Get number of dimensions of input signal */
    nds = mxGetNumberOfDimensions(prhs[0]);
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    if(nds > 2){CANCEL("This method only supports matrix data")}
    /* Check dimensions of weights matrices */
    if ( mxGetM(prhs[1]) != M-1 || mxGetN(prhs[1]) != N )
        {CANCEL("Columns weights matrix must be of size (M-1)xN for an input image of size MxN")}
    if ( mxGetM(prhs[2]) != M || mxGetN(prhs[2]) != N-1 )
        {CANCEL("Columns weights matrix must be of size Mx(N-1) for an input image of size MxN")}

    /* Retrieve input data */
    y = mxGetPr(prhs[0]);
    W1 = mxGetPr(prhs[1]);
    W2 = mxGetPr(prhs[2]);
    
    /* Retrieve rest of options */
    if(nrhs >= 4) ncores = (int)(mxGetPr(prhs[3]))[0];
    else ncores = 1;
    if(nrhs >= 5) maxIters = (int)(mxGetPr(prhs[4]))[0];
    else maxIters = 0;
    
    /* Create output arrays */
    if(nlhs >= 1){
        plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
        x = mxGetPr(plhs[0]);
    }
    else x = (double*)malloc(sizeof(double)*M*N);
    if(nlhs >= 2){
        plhs[1] = mxCreateDoubleMatrix(N_INFO,1,mxREAL);
        info = mxGetPr(plhs[1]);
    }
    else info = NULL;
    
    /* Run solver */
    DR2L1W_TV(M, N, y, W1, W2, x, ncores, maxIters, info);
    
    /* Free resources */
    FREE
    
    return;
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mex.h"
#include "../src/TVopt.h"

/* solveTV2D_PD.cpp

   Solves the 2-dimensial TV proximity problem by applying a Proximal Dykstra method.

   Parameters:
     - 0: multidimensional reference signal y.
     - 1: vector of lambda penalties of each penalty term.
     - 2: vector of dimensions of application of each penalty term.
     - 3: vector of norms of each penalty term.
     - 4: (optional) number of cores to use (default: as defined by environment variable OMP_NUM_THREADS)
     - 5: (optional) maximum number of iterations to run (default: as defined in TVopt.h)
     
   Outputs:
     - 0: solution x.
     - 1: array with optimizer information:
        + [0]: number of iterations run.
        + [1]: stopping tolerance.
*/
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {
    const mwSize *sDims;
    double *x=NULL,*y,*info=NULL,*dims,*norms;
    double *lambdas;
    int *ns=NULL;
    int nds,N,M,npen,ncores,maxIters,i;
    
    #define FREE \
        if(!nlhs) free(x); \
        if(ns) free(ns);
    #define CANCEL(txt) \
        printf("Error in solveTVgen_PDykstrac: %s\n",txt); \
        if(x) free(x); \
        if(info) free(info); \
        return;
        
    /* Check input correctness */
    if(nrhs < 4){CANCEL("not enought inputs");}
    if(!mxIsClass(prhs[0],"double")) {CANCEL("input signal must be in double format")}
    if(!mxIsClass(prhs[1],"double")) {CANCEL("penalties must be in double format")}

    /* Find number of dimensions of the input signal */
    nds=mxGetNumberOfDimensions(prhs[0]);
    /* Get dimensions of input signal */
    sDims = mxGetDimensions(prhs[0]);
    /* Convert dimensions to C array */
    ns = (int*)malloc(sizeof(int)*nds);
    if(!ns) {CANCEL("out of memory")}
    for(i=0;i<nds;i++) ns[i] = (int)sDims[i];
    
    /* Get input signal */
    y = mxGetPr(prhs[0]);
    
    /* Get rest of inputs */
    lambdas = mxGetPr(prhs[1]);
    dims = mxGetPr(prhs[2]);
    norms = mxGetPr(prhs[3]);
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
    npen = (M > N) ? M : N;
    if(nrhs >= 5) ncores = (int)(mxGetPr(prhs[4]))[0];
    else ncores = 1;
    if(nrhs >= 6) maxIters = (int)(mxGetPr(prhs[5]))[0];
    else maxIters = 0;

    /* Create output arrays */
    if(nlhs >= 1){
        plhs[0]=mxCreateNumericArray(nds,sDims,mxDOUBLE_CLASS,mxREAL);
        if(!plhs[0]){CANCEL("out of memory")}
        x=mxGetPr(plhs[0]);
    }
    else x = (double*)malloc(sizeof(double)*mxGetNumberOfElements(prhs[0]));
    if(nlhs >= 2){
        plhs[1] = mxCreateDoubleMatrix(N_INFO,1,mxREAL);
        info = mxGetPr(plhs[1]);
    }
    else info = NULL;
    
    /* Run algorithm */
    PD2_TV(y,lambdas,norms,dims,x,info,ns,nds,npen,ncores,maxIters);
    
    /* Free resources */
    FREE
    
    return;
    
    #undef FREE
    #undef CANCEL
}



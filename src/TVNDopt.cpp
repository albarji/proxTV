/**
    Optimizers for problems dealing with N dimensional versions of TV norm regularization.
    One dimensional solvers are used as building blocks within a proximal stacking framework.
    The proximal combiners are coded here. A function to compute the value of a multidimensional
    TV regularizer is also provided.
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <wchar.h>
#include <time.h>
#include "TVopt.h"
#include "TVmacros.h"

#define LAPACK_ILP64

/*  PD_TV

    Given a reference multidimensional signal y and a series of penalty terms P(x,lambda,d,p), solves the generalized Total Variation
    proximity operator
    
        min_x 0.5 ||x-y||^2 + sum_i P(x,lambda_i,d_i,p_i) .
        
    where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional slice of x following the dimension d_i, TV(x,p) the TV-Lp prox operator for x.
    
    This general formulation allows to apply different TV regularizers over each dimension of the signal in a modular way.
        
    To solve the problem, a Parallel Dykstra-like splitting method is applied, using as building blocks the 1-dimensional TV solvers.
    
    Inputs:
        - y: reference multidimensional signal in vectorized form.
        - lambdas: array of penalty parameters.
        - norms: array indicating the norm for each penalty term.
        - dims: array indicating the dimension over which each penalty term is applied (1 to N).
        - x: array in which to store the solution in vectorized form.
        - info: array in which to store optimizer information (or NULL if not required)
        - ns: array with the lengths of each dimension of y (and x).
        - nds: number of dimensions.
        - npen: number of penalty terms.
        - ncores: number of cores to use
*/
int PD_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters){
    long n,k,nBytes,idx1,idx2;
    int i,j,d,maxDim,iters,nThreads;
    double stop;
    double *xLast=NULL,*out1d=NULL;
    double **p=NULL,**z=NULL;
    long *incs=NULL,*nSlices=NULL;
    Workspace **ws=NULL;

    #define FREE \
        if(p) { for(i=0;i<npen;i++) free(p[i]); free(p); } \
        if(z) { for(i=0;i<npen;i++) free(z[i]); free(z); } \
        if(xLast) free(xLast); \
        if(incs) free(incs); \
        if(nSlices) free(nSlices); \
        if(out1d) free(out1d); \
        if(ws) freeWorkspaces(ws,nThreads);
        
    #define CANCEL(txt,info) \
        printf("PD_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Set number of threads */
    nThreads = (ncores > 1) ? ncores : 1;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Threads: %d\n",nThreads);
    #endif
    omp_set_num_threads(nThreads);
    
    /* Set number of iterations */
    if(maxIters <= 0) maxIters = MAX_ITERS_PD;

    /* Compute total number of elements and maximum along dimensions */
    n = 1; maxDim = 0;
    for(i=0;i<nds;i++){
        n *= ns[i];
        if(ns[i] > maxDim) maxDim = ns[i];
    }
    nBytes = sizeof(double)*n;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"ns: ");
        for(i=0;i<nds;i++)
            fprintf(DEBUG_FILE,"%d ",ns[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"n=%ld\n",n);
        fprintf(DEBUG_FILE,"nBytes=%ld, %d\n",nBytes,(size_t)nBytes);
    #endif
    
    /* Scale penalties */
    for(i=0;i<npen;i++)
        lambdas[i] *= npen;
    
    /* Alloc auxiliary memory */
    p = (double**)calloc(npen,sizeof(double*));
    z = (double**)calloc(npen,sizeof(double*));
    if(!p || !z) {CANCEL("out of memory",info)}
    for(i=0;i<npen;i++){
        p[i] = (double*)malloc(nBytes);
        z[i] = (double*)malloc(nBytes);
        if(!p[i] || !z[i]) {CANCEL("out of memory",info)}
    }
    xLast = (double*)malloc(nBytes);
    incs = (long*)malloc(sizeof(long)*nds);
    nSlices = (long*)malloc(sizeof(long)*nds);
    ws = newWorkspaces(maxDim,nThreads);
    if(!xLast || !incs || !nSlices || !ws) {CANCEL("out of memory",info)}
    
    #ifdef DEBUG
        for(i=0;i<n;i++)
            fprintf(DEBUG_FILE,"%lf ",y[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
    
    /* Initialization */
    #pragma omp parallel for shared(n,x,npen,z,y) private(i,j) default(none)  
    for(i=0;i<n;i++){
        x[i] = 0;
        for(j=0;j<npen;j++)
            z[j][i] = y[i];
    }
        
    /* Computes increments and number of slices to slice the signal over each dimension */
    incs[0] = 1;
    nSlices[0] = n / ns[0];
    for(i=1;i<nds;i++){
        incs[i] = incs[i-1]*ns[i-1];
        nSlices[i] = n / ns[i];
    }
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"incs: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",incs[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"nSlices: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",nSlices[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
    
    /* Main loop */
    stop = DBL_MAX; iters = 0;
    while(stop > STOP_PD && iters < maxIters){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"----------Iteration %d start----------\n",iters+1);
        #endif

        /* Copy actual solution */
        #pragma omp parallel for shared(x,xLast,n) private(i) default(none)  
        for(i=0;i<n;i++){
            xLast[i] = x[i];
            x[i] = 0;
        }
        
        /* Parallelize */
        #pragma omp parallel shared(ws,nSlices,ns,incs,x,p,lambdas,z,norms,npen,dims,DEBUG_FILE) private(d,i,j,k,idx1,idx2) default(none)  
        {
            /* Get thread number */
            int id = omp_get_thread_num();
            Workspace *wsi = ws[id];
            
            /* Prox step for every penalty term */    
            for(i=0;i<npen;i++){
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"··········Penalty %d··········\n",i);
                #endif
                d = dims[i]-1;
                
                int top=nSlices[d];
                wsi->warm = 0;
                
                /* Run 1-dimensional prox operator over each 1-dimensional slice along the specified dimension */    
                #pragma omp for nowait
                for(j=0;j<top;j++){
                    /* Find slice starting point */
                    idx1 = (j / incs[d])*incs[d]*ns[d] + (j % incs[d]);
            
                    /* Construct slice */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        wsi->in[k] = z[i][idx1+idx2];
                    
                    #ifdef DEBUG
                    {
                        int dbgi;
                        fprintf(DEBUG_FILE,"Slice %d: ",j);
                        for(dbgi=0;dbgi<ns[d];dbgi++)
                            fprintf(DEBUG_FILE,"%lf ",wsi->in[dbgi]);
                        fprintf(DEBUG_FILE,"\n");
                    }
                    #endif
                
                    /* Apply 1-dimensional solver */
                    resetWorkspace(wsi);
                    TV(wsi->in, lambdas[i], wsi->out, NULL, ns[d], norms[i], wsi);
                
                    /* Plug solution back */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        p[i][idx1+idx2] = wsi->out[k];
                }
            }
        }
        
        /* Reconstruct signal x from partial penalties solutions */
        for(k=0;k<n;k++)
            for(i=0;i<npen;i++)        
                x[k] += p[i][k] / npen;
    
        /* Z update step */
        #pragma omp parallel for shared(n,npen,z,x,p) private(k,i) default(none)  
        for(k=0;k<n;k++)
            for(i=0;i<npen;i++)
                z[i][k] += x[k] - p[i][k];
        
        /* Compute stopping criterion: mean change */
        stop = 0;
        #pragma omp parallel for shared(x,xLast,n) private(k) reduction(+:stop) default(none)  
        for(k=0;k<n;k++)
            stop += fabs(x[k]-xLast[k]);
        stop /= n;
        
        iters++;
    }
    
    /* Gather output information */
    if(info){
        info[INFO_ITERS] = iters;
        info[INFO_GAP] = stop;
    }
    
    /* Termination check */
    if(iters >= MAX_ITERS_PD){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(PD_TV) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_PD);
        #endif
        if(info) info[INFO_RC] = RC_ITERS;
    }
    else if(info) info[INFO_RC] = RC_OK;
    
    FREE
    return 1;
    
    #undef FREE
    #undef CANCEL
}

/*  PDR_TV

    Given a reference multidimensional signal y and a series of penalty terms P(x,lambda,d,p), solves the generalized Total Variation
    proximity operator
    
        min_x 0.5 ||x-y||^2 + sum_i P(x,lambda_i,d_i,p_i) .
        
    where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional slice of x following the dimension d_i, TV(x,p) the TV-Lp prox operator for x.
    
    This general formulation allows to apply different TV regularizers over each dimension of the signal in a modular way.
        
    To solve the problem, a Parallel Douglas-Rachford splitting method is applied, using as building blocks the 1-dimensional TV solvers.
    
    Inputs:
        - y: reference multidimensional signal in vectorized form.
        - lambdas: array of penalty parameters.
        - norms: array indicating the norm for each penalty term.
        - dims: array indicating the dimension over which each penalty term is applied (1 to N).
        - x: array in which to store the solution in vectorized form.
        - info: array in which to store optimizer information (or NULL if not required)
        - ns: array with the lengths of each dimension of y (and x).
        - nds: number of dimensions.
        - npen: number of penalty terms.
        - ncores: number of cores to use
*/

int PDR_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters){
    long n,k,nBytes,idx1,idx2;
    int i,j,d,maxDim,iters,nThreads;
    double stop;
    double* q;
    double *xLast=NULL,*out1d=NULL;
    double **p=NULL,**z=NULL;
    long *incs=NULL,*nSlices=NULL;
    Workspace **ws=NULL;

    #define FREE \
        if(p) { for(i=0;i<npen;i++) free(p[i]); free(p); } \
        if(q) { free(q); } \
        if(z) { for(i=0;i<npen;i++) free(z[i]); free(z); } \
        if(xLast) free(xLast); \
        if(incs) free(incs); \
        if(nSlices) free(nSlices); \
        if(out1d) free(out1d); \
        if(ws) freeWorkspaces(ws,nThreads);
        
    #define CANCEL(txt,info) \
        printf("PDR_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Set number of threads */
    nThreads = (ncores > 1) ? ncores : 1;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Threads: %d\n",nThreads);
    #endif
    omp_set_num_threads(nThreads);
    
    /* Set number of iterations */
    if(maxIters <= 0) maxIters = MAX_ITERS_DR;

    /* Compute total number of elements and maximum along dimensions */
    n = 1; maxDim = 0;
    for(i=0;i<nds;i++){
        n *= ns[i];
        if(ns[i] > maxDim) maxDim = ns[i];
    }
    nBytes = sizeof(double)*n;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"ns: ");
        for(i=0;i<nds;i++)
            fprintf(DEBUG_FILE,"%d ",ns[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"n=%ld\n",n);
        fprintf(DEBUG_FILE,"nBytes=%ld, %d\n",nBytes,(size_t)nBytes);
    #endif
    
    /* Scale penalties */
    for(i=0;i<npen;i++)
        lambdas[i] *= npen;
    
    /* Alloc auxiliary memory */
    q = (double*)calloc(nBytes,sizeof(double)); // check for out of memory here too
    p = (double**)calloc(npen,sizeof(double*));
    z = (double**)calloc(npen,sizeof(double*));
    if(!q || !p || !z) {CANCEL("out of memory",info)}
    for(i=0;i<npen;i++){
        p[i] = (double*)malloc(nBytes);
        z[i] = (double*)malloc(nBytes);
        if(!p[i] || !z[i]) {CANCEL("out of memory",info)}
    }
    xLast = (double*)malloc(nBytes);
    incs = (long*)malloc(sizeof(long)*nds);
    nSlices = (long*)malloc(sizeof(long)*nds);
    ws = newWorkspaces(maxDim,nThreads);
    if(!xLast || !incs || !nSlices || !ws) {CANCEL("out of memory",info)}
    
    #ifdef DEBUG
        for(i=0;i<n;i++)
            fprintf(DEBUG_FILE,"%lf ",y[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
    
    /* Initialization */
    #pragma omp parallel for shared(n,x,npen,z,y) private(i,j) default(none)  
    for(i=0;i<n;i++){
      x[i]=y[i]/npen; 
      for(j=0;j<npen;j++)
        z[j][i] = y[i];
    }
        
    /* Computes increments and number of slices to slice the signal over each dimension */
    incs[0] = 1;
    nSlices[0] = n / ns[0];
    for(i=1;i<nds;i++){
        incs[i] = incs[i-1]*ns[i-1];
        nSlices[i] = n / ns[i];
    }
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"incs: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",incs[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"nSlices: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",nSlices[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
    
    /* Main loop */
    iters = 0;
    while(iters < maxIters){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"----------Iteration %d start----------\n",iters+1);
        #endif

        /* Copy actual solution */
        #pragma omp parallel for shared(x,xLast,n,q) private(i) default(none)  
        for(i=0;i<n;i++){
            xLast[i] = x[i];
            x[i] = 0;
            q[i]=0;
        }
        
        /* Parallelize */
#pragma omp parallel shared(ws,nSlices,ns,incs,x,p,q,lambdas,z,norms,npen,dims,DEBUG_FILE) private(d,i,j,k,idx1,idx2) default(none)  
        {
            /* Get thread number */
            int id = omp_get_thread_num();
            Workspace *wsi = ws[id];
            
            /* Prox step for every penalty term */    
            for(i=0;i<npen;i++){
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"··········Penalty %d··········\n",i);
                #endif
                d = dims[i]-1;
                
                int top=nSlices[d];
                wsi->warm = 0;
                
                /* Run 1-dimensional prox operators over each 1-dimensional slice along the specified dimension */    
                #pragma omp for nowait
                for(j=0;j<top;j++){
                    /* Find slice starting point */
                    idx1 = (j / incs[d])*incs[d]*ns[d] + (j % incs[d]);
            
                    /* Construct slice */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        wsi->in[k] = z[i][idx1+idx2];
                    
                    #ifdef DEBUG
                    {
                        int dbgi;
                        fprintf(DEBUG_FILE,"Slice %d: ",j);
                        for(dbgi=0;dbgi<ns[d];dbgi++)
                            fprintf(DEBUG_FILE,"%lf ",wsi->in[dbgi]);
                        fprintf(DEBUG_FILE,"\n");
                    }
                    #endif
                
                    /* Apply 1-dimensional solver */
                    resetWorkspace(wsi);
                    if(norms[i] == 1)
                        TV1D_denoise(wsi->in, wsi->out, ns[d], lambdas[i]);
                    else if(norms[i] == 2)
                        morePG_TV2(wsi->in,lambdas[i],wsi->out,NULL,ns[d],wsi);
                    else // general p norm
                        GP_TVp(wsi->in,lambdas[i],wsi->out,NULL,ns[d],norms[i],wsi);
                
                    /* Plug solution back */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        p[i][idx1+idx2] = wsi->out[k];
                }
            }
        }
        
        /* At this point p[i][k] corrs to y[i][n] in (27.25) of Bauschke 
           px[n] will be average of x's
           qy[n] will be average of the prox values
         */
        for(k=0;k<n;k++) {
          for(i=0;i<npen;i++) {
            q[k] += p[i][k] / npen;
            x[k] += z[i][k] / npen;
          }
        }
    
        /* X update step (X for Bauschke, Z for Alvaro) */
        /* X = X + 2Q - P - Y */
#pragma omp parallel for shared(n,npen,z,x,p,q) private(k,i) default(none)  
        for(k=0;k<n;k++)
          for(i=0;i<npen;i++)
            z[i][k] += 2*q[k] - x[k] - p[i][k];
        
        /* Compute stopping criterion: mean change */
        stop = 0;
        #pragma omp parallel for shared(x,xLast,n) private(k) reduction(+:stop) default(none)  
        for(k=0;k<n;k++)
            stop += fabs(x[k]-xLast[k]);
        stop /= n;
        
        iters++;
    }
    
    /* Gather output information */
    if(info){
        info[INFO_ITERS] = iters;
        info[INFO_GAP] = stop;
    }
    
    /* Termination check */
    if(iters >= MAX_ITERS_DR){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(PDR_TV) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_DR);
        #endif
        if(info) info[INFO_RC] = RC_ITERS;
    }
    else if(info) info[INFO_RC] = RC_OK;
    
    FREE
    return 1;
    
    #undef FREE
    #undef CANCEL
}

/*  TVval

    Given a reference multidimensional signal y and a series of penalty terms P(x,lambda,d,p), computes the value of the generalized
    Total Variation norm
    
        sum_i P(x,lambda_i,d_i,p_i) .
        
    where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional slice of x following the dimension d_i, TV(x,p) the TV-Lp norm for x.
    
    Inputs:
        - x: multidimensional signal in vectorized form.
        - lambdas: array of penalty parameters.
        - norms: array indicating the norm for each penalty term.
        - dims: array indicating the dimension over which each penalty term is applied (1 to N).
        - ns: array with the lengths of each dimension of y (and x).
        - nds: number of dimensions.
        - npen: number of penalty terms.
        - ncores: number of cores to use
        
    Output:
        double with value of the generalized TV norm.
*/
double TVval(double *x,double *lambdas,double *norms,double *dims,int *ns,int nds,int npen,int ncores){
    long n,k,nBytes,idx1,idx2;
    int i,j,d,maxDim,iters,nThreads;
    double totalVal;
    long *incs=NULL,*nSlices=NULL;
    Workspace **ws=NULL;

    #define FREE \
        if(incs) free(incs); \
        if(nSlices) free(nSlices); \
        if(ws) freeWorkspaces(ws,nThreads);
        
    #define CANCEL(txt) \
        printf("TVval: %s\n",txt); \
        FREE \
        return 0;
        
    /* Set number of threads */
    nThreads = (ncores > 1) ? ncores : 1;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Threads: %d\n",nThreads);
    #endif
    omp_set_num_threads(nThreads);

    /* Compute total number of elements and maximum along dimensions */
    n = 1; maxDim = 0;
    for(i=0;i<nds;i++){
        n *= ns[i];
        if(ns[i] > maxDim) maxDim = ns[i];
    }
    nBytes = sizeof(double)*n;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"ns: ");
        for(i=0;i<nds;i++)
            fprintf(DEBUG_FILE,"%d ",ns[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"n=%ld\n",n);
        fprintf(DEBUG_FILE,"nBytes=%ld, %d\n",nBytes,(size_t)nBytes);
        fprintf(DEBUG_FILE,"npen=%d\n",npen);
        for(i=0;i<npen;i++)
            fprintf(DEBUG_FILE,"lambda=%lf along dim %d with p=%lf\n",lambdas[i],(int)dims[i],norms[i]);
    #endif
    
    /* Alloc auxiliary memory */
    incs = (long*)malloc(sizeof(long)*nds);
    nSlices = (long*)malloc(sizeof(long)*nds);
    ws = newWorkspaces(maxDim,nThreads);
    if(!incs || !nSlices || !ws) {CANCEL("out of memory")}
    
    #ifdef DEBUG
        for(i=0;i<n;i++)
            fprintf(DEBUG_FILE,"%lf ",x[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
        
    /* Computes increments and number of slices to slice the signal over each dimension */
    incs[0] = 1;
    nSlices[0] = n / ns[0];
    for(i=1;i<nds;i++){
        incs[i] = incs[i-1]*ns[i-1];
        nSlices[i] = n / ns[i];
    }
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"incs: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",incs[i]);
        fprintf(DEBUG_FILE,"\n");
        fprintf(DEBUG_FILE,"nSlices: ");
        for(i=0;i<nds;i++) fprintf(DEBUG_FILE,"%ld ",nSlices[i]);
        fprintf(DEBUG_FILE,"\n");
    #endif
    
        
    /* Parallelize calculation of value */
    #pragma omp parallel shared(x,ws,nSlices,ns,incs,lambdas,norms,npen,dims,DEBUG_FILE) private(d,i,j,k,idx1,idx2) default(none)  
    {
        /* Get thread number */
        int id = omp_get_thread_num();
        Workspace *wsi = ws[id];
        double val;
        
        /* Initialize total TV value for this thread */
        wsi->d[0][0] = 0;
        
        /* Value for every penalty term */    
        for(i=0;i<npen;i++){
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"··········Penalty %d··········\n",i);
            #endif
            d = dims[i]-1;
            
            int top=nSlices[d];
            
            /* Compute 1-dimensional TV norm over each 1-dimensional slice along the specified dimension */    
            #pragma omp for nowait
            for(j=0;j<top;j++){
                /* Find slice starting point */
                idx1 = (j / incs[d])*incs[d]*ns[d] + (j % incs[d]);
        
                /* Construct slice */
                for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                    wsi->in[k] = x[idx1+idx2];
                    
                /* Compute TV-Lp norm value */
                val = 0;
                for(k=0;k<ns[d]-1;k++)
                    val += pow(fabs(wsi->in[k+1] - wsi->in[k]),norms[i]);
                val = lambdas[i] * pow(val,1/norms[i]);
                
                /* Store in workspace structure */
                wsi->d[0][0] += val;
            }
        }
    }
    
    /* Compute total value */
    totalVal = 0;
    for(i=0;i<nThreads;i++)
        totalVal += ws[i]->d[0][0];
    
    FREE
    return totalVal;
    
    #undef FREE
    #undef CANCEL
}

/*  Yang3_TV

    Application of Yags's et al ADMM method for 3D TV-L1.
    
        min_X 0.5 ||X-Y||^2 + lambda * TV_cols(X) + lambda * TV_rows(X) + lambda * TV_depths(X)
        
    The method is bases on calling a 1D TV-L1 solver, for which the standard one provided in this library is used.
    
    Reference: Yang et al - An Efficient ADMM Algorithm for Multidimensional Anisotropic Total Variation Regularization Problems
    
    ADMM parameters are set as recommended in the reference paper, that is

        rho = 10
    
    Inputs:
        int M -- num of rows in input tensor Y
        int N -- num of cols in input tensor Y
        int O -- num of "depths" (third dimension) in input tensor Y
        double* Y -- input tensor
        double lambda    -- regularization weight
        int maxit -- maximum number of iterations to run
        
    Outputs:
        X -- array containing primal solution (pre-alloced space assumed)
        info -- information structure about optimization
*/
int Yang3_TV(size_t M, size_t N, size_t O, double*Y, double lambda, double*X, int maxit, double* info) {
    double *U1 = NULL, *U2 = NULL, *U3 = NULL, *Z1 = NULL, *Z2 = NULL, *Z3 = NULL;
    double rho;
    int i, j, k, idx;
    Workspace *ws = NULL;
    
    #define FREE \
        if(U1) free(U1); \
        if(U2) free(U2); \
        if(U3) free(U3); \
        if(Z1) free(Z1); \
        if(Z2) free(Z2); \
        if(Z3) free(Z3); \
        if(ws) freeWorkspace(ws);
    
    #define CANCEL(txt,info) \
        printf("Yang3_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
    
    // Compute values for optimization parameters
    rho = 10;
    
    // Alloc memory
    int size = (M > N) ? M : N; size = (O > size) ? O : size;
    int totalSize = M*N*O;
    U1 = (double*)calloc(totalSize,sizeof(double));
    U2 = (double*)calloc(totalSize,sizeof(double));
    U3 = (double*)calloc(totalSize,sizeof(double));
    Z1 = (double*)malloc(sizeof(double)*totalSize);
    Z2 = (double*)malloc(sizeof(double)*totalSize);
    Z3 = (double*)malloc(sizeof(double)*totalSize);
    ws = newWorkspace(size);
    if ( !U1 || !U2 || !U3 || !Z1 || !Z2 || !Z3 || !ws )
        {CANCEL("insufficient memory",info)}
        
    // Initialization
    memcpy(Z1, Y, sizeof(double)*totalSize);
    memcpy(Z2, Y, sizeof(double)*totalSize);
    memcpy(Z3, Y, sizeof(double)*totalSize);
    memcpy(X, Y, sizeof(double)*totalSize);
    
    // Set iterations
    if ( maxit <= 0 )
        maxit = MAX_ITERS_YANG;
    
    // Main loop
    int it;
    for ( it = 1 ; it <= maxit ; it++ ) {
        // Update X
        for ( i = 0 ; i < totalSize ; i++ )
            X[i] = (Y[i] + U1[i] + U2[i] + U3[i] + rho*(Z1[i] + Z2[i] + Z3[i])) / (1 + 3*rho);
            
        // Update Z1: prox operators along first dimension (column fibers)
        for ( i = 0 ; i < N ; i++ ) {
            for ( j = 0 ; j < O ; j++ ) {
                // Copy first dimension fiber to workspace
                for ( k = 0 ; k < M ; k++ ) {
                    idx = k + M * ( i + N * j );
                    ws->in[k] = -1. / rho * U1[idx] + X[idx];
                }
                resetWorkspace(ws);
                TV(ws->in, lambda/rho, ws->out, NULL, M, 1, ws);
                // Recover data
                memcpy(Z1 + M * ( i + N * j ), ws->out, sizeof(double)*M);
            }
        }
        
        // Update Z2: prox operators along second dimension (row fibers)
        for ( k = 0 ; k < M ; k++ ) {
            for ( j = 0 ; j < O ; j++ ) {
                // Copy second dimension fiber to workspace
                for ( i = 0 ; i < N ; i++ ) {
                    idx = k + M * ( i + N * j );
                    ws->in[i] = -1. / rho * U2[idx] + X[idx];
                }
                resetWorkspace(ws);
                TV(ws->in, lambda/rho, ws->out, NULL, N, 1, ws);
                // Recover data
                for ( i = 0 ; i < N ; i++ ) {
                    idx = k + M * ( i + N * j );
                    Z2[idx] = ws->out[i];
                }
            }
        }
        
        // Update Z3: prox operators along third dimension (depth fibers)
        for ( k = 0 ; k < M ; k++ ) {
            for ( i = 0 ; i < N ; i++ ) {
                // Copy third dimension fiber to workspace
                for ( j = 0 ; j < O ; j++ ) {
                    idx = k + M * ( i + N * j );
                    ws->in[j] = -1. / rho * U3[idx] + X[idx];
                }
                resetWorkspace(ws);
                TV(ws->in, lambda/rho, ws->out, NULL, O, 1, ws);
                // Recover data
                for ( j = 0 ; j < O ; j++ ) {
                    idx = k + M * ( i + N * j );
                    Z3[idx] = ws->out[j];
                }
            }
        }
        
        // Update Us
        for ( i = 0 ; i < totalSize ; i++ ) {
            U1[i] += rho * (Z1[i] - X[i]);
            U2[i] += rho * (Z2[i] - X[i]);
            U3[i] += rho * (Z3[i] - X[i]);
        }
    }
    
    // Gather output information
    if(info){
        info[INFO_ITERS] = it;
        info[INFO_RC] = RC_OK;
    }
    
    // Free memory and return
    FREE
    
    #undef FREE
    #undef CANCEL
}


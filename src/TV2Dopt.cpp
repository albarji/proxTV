/**
    Optimizers for problems dealing with 2 dimensional versions of TV norm regularization.
    One dimensional solvers are used as building blocks within a proximal stacking framework.
    The proximal combiners are coded here.

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

/* Internal functions */
void DR_proxDiff(size_t n, double* input, double* output, double W, double norm, Workspace *ws);
void DR_columnsPass(size_t M, size_t N, double* input, double* output, double W, double norm, Workspace **ws);
void DR_rowsPass(size_t M, size_t N, double* input, double* output, double* ref, double W, double norm, Workspace **ws);

/* Internal definitions */
#define ALG_CONDAT 0
#define ALG_CHAMBOLLE_POCK 1
#define ALG_CHAMBOLLE_POCK_ACC 2

/*  PD2_TV

    Optimized version of Proximal Dykstra for the case where only 1 or 2 penalty terms appear.

    Given a reference multidimensional signal y and a series of penalty terms P(x,lambda,d,p), solves the generalized Total Variation
    proximity operator

        min_x 0.5 ||x-y||^2 + sum_{i=1}^2 P(x,lambda_i,d_i,p_i) .

    where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional slice of x following the dimension d_i, TV(x,p) the TV-Lp prox operator for x.

    This general formulation allows to apply different TV regularizers over each dimension of the signal in a modular way.

    To solve the problem, a Dykstra-like splitting method is applied, using as building blocks the 1-dimensional TV solvers.

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
int PD2_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters){
    long n,k,nBytes,idx1,idx2;
    int i,j,d,maxDim,iters;
    double stop;
    double *xLast=NULL;
    double *p=NULL,*z=NULL,*q=NULL;
    long *incs=NULL,*nSlices=NULL;
    Workspace **ws=NULL;
    short nThreads;

    #define FREE \
        if(p) free(p); \
        if(z) free(z); \
        if(q) free(q); \
        if(xLast) free(xLast); \
        if(incs) free(incs); \
        if(nSlices) free(nSlices); \
        if(ws) freeWorkspaces(ws,nThreads);

    #define CANCEL(txt,info) \
        printf("PD2_TV: %s\n",txt); \
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

    /* This algorithm can only work with 1 or 2 penalties */
    if(npen > 2)
        {CANCEL("this algorithm can not work with more than 2 penalties",info)}

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

    /* Alloc auxiliary memory */
    p = (double*)malloc(nBytes);
    z = (double*)malloc(nBytes);
    q = (double*)malloc(nBytes);
    if(!p || !z || !q) {CANCEL("out of memory",info)}
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
    #pragma omp parallel for shared(x,y,p,q,n) private(i) default(none)
    for(i=0;i<n;i++){
        x[i] = y[i];
        p[i] = 0;
        q[i] = 0;
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
    while(stop > STOP_PD && (npen > 1 || !iters) && iters < maxIters){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"----------Iteration %d start----------\n",iters+1);
        #endif

        /* Copy actual solution */
        #pragma omp parallel for shared(x,xLast,n) private(i) default(none)
        for(i=0;i<n;i++)
            xLast[i] = x[i];

        /* Prox step for the first penalty term */
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"··········Penalty 0··········\n");
        #endif
        d = int(dims[0]-1);
        /* Run 1-dimensional prox operator over each 1-dimensional slice along the specified dimension (parallelized) */
        #pragma omp parallel shared(ws,nSlices,ns,d,incs,x,p,lambdas,z,norms) private(j,k,idx1,idx2) default(none)
        {
            /* Get thread number */
            int id = omp_get_thread_num();
            Workspace *wsi = ws[id];
            wsi->warm = 0;
            int top=nSlices[d];

            #pragma omp for
            for(j=0;j<top;j++){
                /* Find slice starting point */
                idx1 = (j / incs[d])*incs[d]*ns[d] + (j % incs[d]);

                /* Construct slice */
                for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                    wsi->in[k] = x[idx1+idx2]+p[idx1+idx2];

                /* Apply 1-dimensional solver */
                resetWorkspace(wsi);
                TV(wsi->in, lambdas[0], wsi->out, NULL, ns[d], norms[0], wsi);

                /* Plug solution back */
                for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                    z[idx1+idx2] = wsi->out[k];
            }
        }

        /* Update p */
        #pragma omp parallel for shared(p,x,z,n) private(i) default(none)
        for(i=0;i<n;i++)
            p[i] += x[i] - z[i];

        /* Prox step for the second penalty term (if any) */
        if(npen >= 2){
            d = int(dims[1]-1);

            /* Run 1-dimensional prox operator over each 1-dimensional slice along the specified dimension (parallelized) */
            #pragma omp parallel shared(ws,nSlices,ns,d,incs,x,q,lambdas,z,norms) private(j,k,idx1,idx2) default(none)
            {
                /* Get thread number */
                int id = omp_get_thread_num();
                Workspace *wsi = ws[id];
                wsi->warm = 0;
                int top=nSlices[d];

                #pragma omp for
                for(j=0;j<top;j++){
                    /* Find slice starting point */
                    idx1 = (j / incs[d])*incs[d]*ns[d] + (j % incs[d]);

                    /* Construct slice */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        wsi->in[k] = z[idx1+idx2] + q[idx1+idx2];

                    /* Apply 1-dimensional solver */
                    resetWorkspace(wsi);
                    TV(wsi->in, lambdas[1], wsi->out, NULL, ns[d], norms[1], wsi);

                    /* Plug solution back */
                    for(k=0,idx2=0 ; k<ns[d] ; k++,idx2+=incs[d])
                        x[idx1+idx2] = wsi->out[k];
                }
            }

            /* Update q */
            #pragma omp parallel for shared(q,z,x,n) private(i) default(none)
            for(i=0;i<n;i++)
                q[i] += z[i] - x[i];
        }
        else{
            #pragma omp parallel for shared(z,x,n) private(i) default(none)
            for(i=0;i<n;i++)
                x[i] = z[i];
            memcpy(x,z,nBytes);
        }

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
        if(info) info[INFO_RC] = RC_ITERS;
    }
    else if(info) info[INFO_RC] = RC_OK;

    FREE
    return 1;

    #undef FREE
    #undef CANCEL
}

/* DR2_TV

    Applies a variant of the Douglas Rachford algorithm, a.k.a. Alternating
    Reflections method, to solve the 2D TV proximity problem.

    To do so the 2D TV proximity problem

        min_x 0.5 ||x-y||^2 + TV_cols(x) + TV_rows(x)

    is rewritten as the equivalent problem.

        min_x 0.5 ||x||^2 + TV_cols(x) + TV_rows*(x)

    where TV_rows*(x) = TV_rows(x) + x^T y. This problem takes the
    form of the minimization of two Lovasz extensions, which in turn
    can be shown to be solvable by the following best approximation problem:

        min_{y1,y2} ||y_1 - y_2|| s.t. y_1 \in B_{-rows*}, y_2 \in B_{cols}

    where B are the base polytopes of the corresponding Lovasz extensions.
    Projection onto these polytopes can be computed easily through prox
    operators of 1D TV, as it can be shown that

        proj[B_{cols}](z) = z - prox_cols(z)
        proj[B_{-rows*}](y) = z + prox_row(-z+y)

    The method thus alternates reflections based on these projections to
    solve the problem.

     Inputs
         int M -- num of rows in input matrix (unary)
         int N -- num of cols in input matrix
         double* unary -- input matrix / image
         double W1    -- regularization weight along cols
         double W2    -- regularization weight along rows
         double norm1 -- TV norm to use along cols
         double norm2 -- TV norm to use along rows
         int nThreads  -- num of threads to use for parallel processing

     Outputs:
         s -- array containing primal solution (pre-alloced space assumed)

    References
        - Jegelka, Bach, Sra - Reflection methods for user-friendly
            submodular optimization
        - Bach - Learning with Submodular Functions: A Convex Optimization
            Perspective
 */
int DR2_TV(size_t M, size_t N, double*unary, double W1, double W2,
  double norm1, double norm2, double*s, int nThreads, int maxit, double* info)
{

  size_t i;
  double *t = NULL;
  double *tb = NULL;
  Workspace **ws = NULL;

  #define FREE \
    if(t) free(t); \
    if(tb) free(tb); \
    if(ws) freeWorkspaces(ws,nThreads);

  #define CANCEL(txt,info) \
      printf("DR2_TV: %s\n",txt); \
      FREE \
      if(info) info[INFO_RC] = RC_ERROR;\
      return 0;

  if (nThreads < 1) nThreads = 1;
  omp_set_num_threads(nThreads);
  int maxDim = int( (M > N) ? M : N );

  // Alloc memory for algorithm */
  t = (double*) malloc(sizeof(double)*M*N);
  tb = (double*)malloc(sizeof(double)*M*N);
  ws = newWorkspaces(maxDim,nThreads);

  if (!t || !tb || !ws)
    {CANCEL("out of memory", info)}

  /* Set number of iterations */
  if(maxit <= 0) maxit = MAX_ITERS_DR;

  // Normalize input to start DR: t = ones(size(unary)) * mean(unary(:));
  double sum=0;
  for (i=0; i < M*N; i++)
    sum += unary[i];
  sum = 2*sum / (M*N);
  for (i=0; i < M*N; i++)
    t[i]=sum;

    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Starting Douglas-Rachford with size=[%d,%d], penalties=[%lf,%lf], norms=[%lf,%lf], threads=%d\n",M,N,W1,W2,norm1,norm2,nThreads); fflush(DEBUG_FILE);
    #endif

  int iter = 0;
  /* MAIN LOOP */
  while ( iter < maxit ) {
    iter++;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Reflection along columns\n"); fflush(DEBUG_FILE);
    #endif
    // Projection for B_{cols}
    DR_columnsPass(M, N, t, s, W1, norm1, ws);
    // Reflection for B_{cols}
    for (i=0; i < M*N; i++) s[i] = 2*s[i] - t[i];

    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Dual projection along rows\n"); fflush(DEBUG_FILE);
    #endif
    // Projection for B_{-rows*}
    DR_rowsPass(M, N, s, tb, unary, W2, norm2, ws);
    // Reflection for B_{-rows*}
    for (i=0; i < M*N; i++) tb[i] = 2*tb[i] - s[i];

    // Combiner step to get t_{i+1} = 0.5 [R[B_{-rows*}] R[B_{cols}] + I] t_i
    for (i=0; i < M*N; i++) t[i] = 0.5*(t[i]+tb[i]);
  }

  // DR is divergent, but with an additional projection we can recover valid solutions
  // Projection for B_{cols}
  DR_columnsPass(M, N, t, s, W1, norm1, ws);
  // Projection for B_{-rows*}
  DR_rowsPass(M, N, s, tb, unary, W2, norm2, ws);
  for (i = 0; i < M*N; i++) s[i] = tb[i] - s[i];

    /* Gather output information */
    if(info){
        info[INFO_ITERS] = iter;
        info[INFO_RC] = RC_OK;
    }

  // Free and return
  FREE
  return 0;

  #undef FREE
  #undef CANCEL
}

/**
    Projection agains the base polyhedron of the columns in the DR algorithm.

    proj[B_{cols}](z) = z - prox_cols(z)

    @param M number of rows in input signal
    @param N number of columns in input signal
    @param input signal
    @param output signal
    @param W weight of TV regularization
    @param norm degree of TV
    @param ws array of Workspaces to use for the computation
*/
void DR_columnsPass(size_t M, size_t N, double* input, double* output, double W, double norm, Workspace **ws) {
    #pragma omp parallel shared(M,N,input,output,W,norm,ws) default(none)
    {
        // Get thread number
        int id = omp_get_thread_num();
        // Get corresponding workspace
        Workspace *wsi = ws[id];
        wsi->warm = 0;

        // Run 1-d solvers in parallel on each column of the input
        #pragma omp for
        for (int j=0; j < N; j++) {
            resetWorkspace(wsi);
            // Prepare inputs
            memcpy(wsi->in, input+(M*j), sizeof(double)*M);
            // Compute prox difference for this column
            DR_proxDiff(M, wsi->in, wsi->out, W, norm, wsi);
            // Save output
            memcpy(output+(M*j), wsi->out, sizeof(double)*M);
        }
    }
}

/**
    Projection agains the base polyhedron of the (modified) rows in the
    DR algorithm.

    proj[B_{-rows*}](z) = z + prox_row(-z+y)

    @param M number of rows in input signal
    @param N number of columns in input signal
    @param input signal
    @param output signal
    @param ref reference signal (input to the global DR algorithm)
    @param W weight of TV regularization
    @param norm degree of TV
    @param ws array of Workspaces to use for the computation
*/
void DR_rowsPass(size_t M, size_t N, double* input, double* output, double* ref, double W, double norm, Workspace **ws) {
    #pragma omp parallel shared(M,N,input,ref,output,W,norm,ws) default(none)
    {
        int i,j;
        // Get thread number
        int id = omp_get_thread_num();
        // Get corresponding workspace
        Workspace *wsi = ws[id];
        wsi->warm = 0;

        // Run 1-d solvers in parallel on each row of the input
        #pragma omp for
        for (j=0; j < M; j++) {
            resetWorkspace(wsi);
            // Prepare inputs, considering displacement from reference signal
            size_t idx;
            for ( idx = j, i = 0 ; i < N ; i++, idx+=M )
                wsi->in[i] = ref[idx]-input[idx];
            // Compute prox difference for this row
            DR_proxDiff(N, wsi->in, wsi->out, W, norm, wsi);
            // Save output, recovering displacement from reference signal
            for ( idx = j, i = 0 ; i < N ; i++, idx+=M )
                output[idx] = ref[idx] - wsi->out[i];
        }
    }
}

/**
 Applies the TV proximity operator to the given 1D input and returns the
difference between input and output of prox.

 This is equivalent to computing projection against the base polytope of
 the corresponding TV function.

 @param n length of input signal
 @param input signal
 @param output signal
 @param W weight of TV regularization
 @param norm degree of TV
 @param ws Workspace to use for the computation
 */
void DR_proxDiff(size_t n, double* input, double* output, double W, double norm, Workspace *) {

    // Compute proximity
    TV(input, W, output, NULL, (int) n , norm, NULL);
    // Return differences between input and proximity output
    for (size_t i=0; i < n; i++)
      output[i] = input[i] - output[i];
}

/*  CondatChambollePock2_TV

    Application of Condat's or Chambolle-Pock generic proximal algorithm for 2D TV-L1.

        min_X 0.5 ||X-Y||^2 + lambda * TV_cols(X) + lambda * TV_rows(X)

    Although both algorithms differ in nature, when applied to 2D TV-L1 they result in almost the same iterations except
    for the step dealing with the differentiable term (||X-Y||^2) which in Chambolle-Pock is processed through a prox operator
    and in Condat's the gradient is used. This results in the different steps:

        Condat: Xt = X - TAU * (X - Y) - TAU * ( ([zeros(1,N);U1] - [U1;zeros(1,N)]) + ([zeros(M,1),U2] - [U2,zeros(M,1)] ));
        Chamboll-Pock: Xt = 1/(1+TAU) * (X + TAU * Y - TAU * ( ([zeros(1,N);U1] - [U1;zeros(1,N)]) + ([zeros(M,1),U2] - [U2,zeros(M,1)] )));

    Algorithm's optimization parameters are tuned in accordance to the values used
    in Toolbox Sparse Optmization by Gabriel Peyre, where an implementation of
    Chambolle-Pock based on ADMM is provided
    (http://www.mathworks.com/matlabcentral fileexchange/16204-toolbox-sparse-optmization).
    The rules for choosing the values seem to make use of the fact that the norm of the
    linear operator in 2D-TV is 8, and so:

        theta = 1  (makes the algorithm analogous to ADMM)
        sigma = 10;    (no explanation for this, but works well)
        tau = 0.9 / (8*sigma);  (the 0.9 seems to come from the fact that we must have tau < 1, and the rest from tau*sigma*L^2 < 1)

    Inputs:
        int M -- num of rows in input matrix Y
        int N -- num of cols in input matrix Y
        double* Y -- input matrix / image
        double lambda    -- regularization weight
        short alg -- algorithm choice
            0: Condat
            1: Chambolle-Pock
        int maxit -- maximum number of iterations to run

    Outputs:
        X -- array containing primal solution (pre-alloced space assumed)
        info -- information structure about optimization
*/
int CondatChambollePock2_TV(size_t M, size_t N, double*Y, double lambda, double*X, short alg, int maxit, double* info) {
    double sigma, tau, theta, normalizer, tmp;
    double gamma = -1;
    double *U1 = NULL, *U2 = NULL, *Xt = NULL, *Z = NULL;
    size_t i, j;

    #define FREE \
        if(U1) free(U1); \
        if(U2) free(U2); \
        if(Xt) free(Xt); \
        if(Z) free(Z);

    #define CANCEL(txt,info) \
        printf("Condat2_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;

    // Check algorithm choice
    if ( alg != ALG_CONDAT && alg != ALG_CHAMBOLLE_POCK && alg != ALG_CHAMBOLLE_POCK_ACC)
        {CANCEL("Algorithm parameter has an invalid value",info)}

    // Compute values for optimization parameters
    if ( alg == ALG_CHAMBOLLE_POCK_ACC ) {
        sigma = 10;
        tau = .9/(sigma*8);
        gamma = 1./lambda;
        theta = 1.;
    } else {
        sigma = 10;
        tau = .9/(sigma*8);
        theta = 1.;
    }

    // Alloc memory
    U1 = (double*)malloc(sizeof(double)*(M-1)*N);
    U2 = (double*)malloc(sizeof(double)*M*(N-1));
    Xt = (double*)malloc(sizeof(double)*M*N);
    Z = (double*)malloc(sizeof(double)*M*N);
    if ( !U1 || !U2 || !Xt || !Z )
        {CANCEL("insufficient memory",info)}

    // Initialize values

    // X = Y
    memcpy(X, Y, sizeof(double)*M*N);
    // U1 = Y(2:end,:) - Y(1:end-1,:);
    for ( i = 0 ; i < M-1 ; i++ ) {
        for ( j = 0 ; j < N ; j++ )
            U1[i+(M-1)*j] = Y[i+1+M*j] - Y[i+M*j];
    }
    // U2 = Y(:,2:end) - Y(:,1:end-1);
    for ( i = 0 ; i < M ; i++ ) {
        for ( j = 0 ; j < N-1 ; j++ )
            U2[i+M*j] = Y[i+M*(j+1)] - Y[i+M*j];
    }

    // Set iterations
    if ( maxit <= 0 )
        maxit = MAX_ITERS_CONDAT;

    // Main loop
    int it; double stop = DBL_MAX;

    for ( it = 1 ; stop > STOP_CONDAT && it <= maxit ; it++ ) {
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Condat iteration %d\n",it); fflush(DEBUG_FILE);
        #endif

        // Smooth term step

        // [zeros(1,N);U1] - [U1;zeros(1,N)]
        for ( j = 0 ; j < N ; j++ ) // First row
            Xt[M*j] = -U1[(M-1)*j];
        for ( i = 1 ; i < M-1 ; i++ ) { // Intermediate rows
            for ( j = 0 ; j < N ; j++ )
                Xt[i+M*j] = U1[i-1+(M-1)*j] - U1[i+(M-1)*j];
        }
        for ( j = 0 ; j < N ; j++ ) // Last row
            Xt[M-1+M*j] = U1[M-2+(M-1)*j];

        // [zeros(M,1),U2] - [U2,zeros(M,1)]
        for ( i = 0 ; i < M ; i++ ) // First column
            Xt[i] += -U2[i];
        for ( i = 0 ; i < M ; i++ ) { // Intermediate columns
            for ( j = 1 ; j < N-1 ; j++ )
                Xt[i+M*j] += U2[i+M*(j-1)] - U2[i+M*j];
        }
        for ( i = 0 ; i < M ; i++ ) // Last column
            Xt[i+M*(N-1)] += U2[i+M*(N-2)];

        // Condat's method: use gradient
        if ( alg == ALG_CONDAT) {
            // Gradient step: Xt = X - TAU * (X - Y) - TAU * ( ([zeros(1,N);U1] - [U1;zeros(1,N)]) + ([zeros(M,1),U2] - [U2,zeros(M,1)] ));

            // X - TAU * (X - Y) - TAU * (result_so_far)
            for ( i = 0 ; i < M*N ; i++ )
                Xt[i] = X[i] - tau * (X[i] - Y[i] + Xt[i]);

        // Chambolle-Pock's method: use proximal operator
        } else if ( alg == ALG_CHAMBOLLE_POCK ||  alg == ALG_CHAMBOLLE_POCK_ACC ) {
            // Proximal step: Xt = 1/(1+TAU) * (X + TAU * Y - TAU * ( ([zeros(1,N);U1] - [U1;zeros(1,N)]) + ([zeros(M,1),U2] - [U2,zeros(M,1)] )));

            // 1/(1+TAU) * (X + TAU * Y - TAU * (result_so_far))
            tmp = 1./(1.+tau);
            for ( i = 0 ; i < M*N ; i++ )
                Xt[i] = tmp * (X[i] + tau * (Y[i] - Xt[i]));
        }

        // Combiner step: Z = Xt + \theta (Xt - X);
        for ( i = 0 ; i < M*N ; i++ )
            Z[i] = Xt[i] + theta * (Xt[i] - X[i]);
        // Update algorithm parameters (only in accelerated Chambolle-Pock)
        if ( alg == ALG_CHAMBOLLE_POCK_ACC ) {
            tau *= theta;
            sigma /= theta;
            theta = 1. / sqrt(1 + 2 * gamma * tau); // 1/sqrt(1 + 2*gamma*tau)
        }

        // Compute stopping tolerance before copying X
        stop = 0; normalizer = 0;
        for ( i = 0 ; i < M*N ; i++ ) {
            normalizer += X[i] * X[i];
            tmp = Xt[i] - X[i];
            stop += tmp*tmp;
        }
        stop = sqrt(stop/normalizer);

        // Store solution at current step: X = Xt;
        memcpy(X, Xt, sizeof(double) * M * N);

        // Proximal (projection) for rows: linf_proj(U1 + SIGMA * (Z(2:end,:) - Z(1:end-1,:)), lambda);

        // U1 + SIGMA * (Z(2:end,:) - Z(1:end-1,:))
        for ( i = 1 ; i < M ; i++ ) {
            for ( j = 0 ; j < N ; j++ )
                U1[i-1+(M-1)*j] += sigma * (Z[i+M*j] - Z[i-1+M*j]);
        }
        // U1 = linf_proj(U1,lambda)
        for ( i = 0 ; i < N*(M-1) ; i++ ) {
            if ( U1[i] < -lambda )
                U1[i] = -lambda;
            else if ( U1[i] > lambda )
                U1[i] = lambda;
        }

        // Proximal (projection) for columns: linf_proj(U2 + SIGMA * (Z(:,2:end) - Z(:,1:end-1)), lambda);

        // U2 + SIGMA * (Z(:,2:end) - Z(:,1:end-1))
        for ( i = 0 ; i < M ; i++ ) {
            for ( j = 1 ; j < N ; j++ )
                U2[i+M*(j-1)] += sigma * (Z[i+M*j] - Z[i+M*(j-1)]);
        }
        // U2 = linf_proj(U2,lambda)
        for ( i = 0 ; i < M*(N-1) ; i++ ) {
            if ( U2[i] < -lambda )
                U2[i] = -lambda;
            else if ( U2[i] > lambda )
                U2[i] = lambda;
        }
    }

    // Gather output information
    if(info){
        info[INFO_ITERS] = it;
        info[INFO_RC] = RC_OK;
    }

    // Free memory and return
    FREE
    return 1;

    #undef FREE
    #undef CANCEL
}

/*  Yang2_TV

    Application of Yags's et al ADMM method for 2D TV-L1.

        min_X 0.5 ||X-Y||^2 + lambda * TV_cols(X) + lambda * TV_rows(X)

    The method is bases on calling a 1D TV-L1 solver, for which the standard one provided in this library is used.

    Reference: Yang et al - An Efficient ADMM Algorithm for Multidimensional Anisotropic Total Variation Regularization Problems

    ADMM parameters are set as recommended in the reference paper, that is

        rho = 10

    Inputs:
        int M -- num of rows in input matrix Y
        int N -- num of cols in input matrix Y
        double* Y -- input matrix / image
        double lambda    -- regularization weight
        int maxit -- maximum number of iterations to run

    Outputs:
        X -- array containing primal solution (pre-alloced space assumed)
        info -- information structure about optimization
*/
int Yang2_TV(size_t M, size_t N, double*Y, double lambda, double*X, int maxit, double* info) {
    double *U1 = NULL, *U2 = NULL, *Z1 = NULL, *Z2 = NULL;
    double rho;
    size_t i, j;
    Workspace *ws = NULL;

    #define FREE \
        if(U1) free(U1); \
        if(U2) free(U2); \
        if(Z1) free(Z1); \
        if(Z2) free(Z2); \
        if(ws) freeWorkspace(ws);

    #define CANCEL(txt,info) \
        printf("Yang2_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;

    // Compute values for optimization parameters
    rho = 10;

    // Alloc memory
    int size = int( (M > N) ? M : N );
    U1 = (double*)calloc(M*N,sizeof(double));
    U2 = (double*)calloc(M*N,sizeof(double));
    Z1 = (double*)malloc(sizeof(double)*M*N);
    Z2 = (double*)malloc(sizeof(double)*M*N);
    ws = newWorkspace(size);
    if ( !U1 || !U2 || !Z1 || !Z2 || !ws )
        {CANCEL("insufficient memory",info)}

    // Initialization
    memcpy(Z1, Y, sizeof(double)*M*N);
    memcpy(Z2, Y, sizeof(double)*M*N);
    memcpy(X, Y, sizeof(double)*M*N);

    // Set iterations
    if ( maxit <= 0 )
        maxit = MAX_ITERS_YANG;

    // Main loop
    int it;
    for ( it = 1 ; it <= maxit ; it++ ) {
        // Update X
        for ( i = 0 ; i < M*N ; i++ )
            X[i] = (Y[i] + U1[i] + U2[i] + rho*(Z1[i] + Z2[i])) / (1 + 2*rho);

        // Update Z1 (prox operators along rows)
        for ( i = 0 ; i < M ; i++ ) {
            // Copy row data to workspace
            for ( j = 0 ; j < N ; j++ )
                ws->in[j] = -1. / rho * U1[j*M+i] + X[j*M+i];
            resetWorkspace(ws);
            TV(ws->in, lambda/rho, ws->out, NULL, (int) N, 1, ws);
            // Recover data
            for ( j = 0 ; j < N ; j++ )
                Z1[j*M+i] = ws->out[j];
        }

        // Update Z2 (prox operators along columns)
        for ( i = 0 ; i < N ; i++ ) {
            // Copy column data to workspace
            for ( j = 0 ; j < M ; j++ )
                ws->in[j] = -1. / rho * U2[i*M+j] + X[i*M+j];
            TV(ws->in, lambda/rho, ws->out, NULL, (int) M, 1, ws);
            // Recover data
            memcpy(Z2+i*M, ws->out, sizeof(double)*M);
        }

        // Update U1
        for ( i = 0 ; i < M*N ; i++ )
            U1[i] += rho * (Z1[i] - X[i]);
        // Update U2
        for ( i = 0 ; i < M*N ; i++ )
            U2[i] += rho * (Z2[i] - X[i]);
    }

    // Gather output information
    if(info){
        info[INFO_ITERS] = it;
        info[INFO_RC] = RC_OK;
    }

    // Free memory and return
    FREE
    return 1;

    #undef FREE
    #undef CANCEL
}

/*  Kolmogorov2_TV

    Application of Kolmogorov et al's proximal algorithm for 2D TV-L1.

        min_X 0.5 ||X-Y||^2 + TV_cols(X) + TV_rows(X)

    This method uses an ADMM-like splitting to deal with rows and columns separately,
    getting the problem

        min_X max_U X·U + TV_rows(X) + 0.5 ||X-Y||^2 - TV_cols*(U)

    where TV_cols* is the Fenchel dual of TV_cols, which is dealt with using Moreau's
    identity.

    Inputs:
        int M -- num of rows in input matrix Y
        int N -- num of cols in input matrix Y
        double* Y -- input matrix / image
        double lambda    -- regularization weight
        int maxit -- maximum number of iterations to run

    Outputs:
        X -- array containing primal solution (pre-alloced space assumed)
        info -- information structure about optimization

    References:
        Kolmogorov et al - Total variation on a Tree [Page 19]
*/
int Kolmogorov2_TV(size_t M, size_t N, double*Y, double lambda, double*X, int maxit, double* info) {
    double sigma, tau, theta, tmp, normalizer;
    double *U = NULL, *Xprev = NULL, *row=NULL, *rowout=NULL, *TMP=NULL;
    size_t i, j;
    size_t NM = N*M;

    #define FREE \
        if(U) free(U); \
        if(Xprev) free(Xprev); \
        if(row) free(row); \
        if(rowout) free(rowout); \
        if(TMP) free(TMP);

    #define CANCEL(txt,info) \
        printf("Kolmogorov2_TV: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;

    // Initial stepsize values
    theta = 1.;
    tau = 1. / 2.;
    sigma = 1.;

    // Alloc memory
    U = (double*)malloc(sizeof(double)*M*N);
    Xprev = (double*)malloc(sizeof(double)*M*N);
    TMP = (double*)malloc(sizeof(double)*M*N);
    row = (double*)malloc(sizeof(double)*N);
    rowout = (double*)malloc(sizeof(double)*N);
    if ( !U || !Xprev || !row || !rowout || !TMP )
        {CANCEL("insufficient memory",info)}

    // Initialize values

    // X = Xprev = U = Y
    memcpy(X, Y, sizeof(double)*NM);
    memcpy(Xprev, Y, sizeof(double)*NM);
    memcpy(U, Y, sizeof(double)*NM);

    // Set iterations
    if ( maxit <= 0 )
        maxit = MAX_ITERS_KOLMOGOROV;

    // Main loop
    int it; double stop = DBL_MAX;

    for ( it = 1 ; stop > STOP_KOLMOGOROV && it <= maxit ; it++ ) {
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Kolmogorov iteration %d\n",it); fflush(DEBUG_FILE);
        #endif

        // Columns step (dual U update)

        // Input to dual prox solver: z = U + sigma*(X + theta *(X-Xprev)); input = sigma^-1*z
        for ( i = 0 ; i < NM ; i++ ) {
            TMP[i] = U[i] + sigma * (X[i] + theta * (X[i] - Xprev[i]));
            TMP[i] /= sigma;
        }
        // Dual prox for cols
        for ( j = 0 ; j < NM ; j+=M ) {
            // Normal prox
            TV(TMP+j, lambda/sigma, U+j, NULL, (int) M, 1, NULL);
            // Moreau's identity
            for ( i = 0 ; i < M ; i++ )
                (U+j)[i] = sigma*((TMP+j)[i] - (U+j)[i]);
        }

        // Save current X
        memcpy(Xprev, X, sizeof(double)*NM);

        // Rows step (primal X update)

        // Input to prox solver: z = X - tau * U; input = (1-tau^-1)^-1*(y+tau^-1*z)
        for ( i = 0 ; i < NM ; i++ ) {
            TMP[i] = X[i] - tau * U[i];
            TMP[i] = 1/(1+1/tau) * (Y[i] + 1/tau * TMP[i]);
        }
        // Prox operator for rows
        for ( i = 0 ; i < M ; i++ ) {
            // Copy input
            for ( j = 0 ; j < N ; j++ )
                row[j] = TMP[j*M+i];
            // Prox operator
            TV(row, lambda/(1.+1./tau), rowout, NULL, (int) N, 1, NULL);
            // Recover output
            for ( j = 0 ; j < N ; j++ )
                X[j*M+i] = rowout[j];
        }

        // Update stepsizes
        theta = 1. / sqrt(1 + 1*tau);
        tau *= theta;
        sigma /= theta;

        // Compute stopping tolerance
        stop = 0; normalizer = 0;
        for ( i = 0 ; i < NM ; i++ ) {
            normalizer += X[i] * X[i];
            tmp = Xprev[i] - X[i];
            stop += tmp*tmp;
        }
        stop = sqrt(stop/normalizer);
    }

    // Gather output information
    if(info){
        info[INFO_ITERS] = it;
        info[INFO_RC] = RC_OK;
    }

    // Free memory and return
    FREE
    return 1;

    #undef FREE
    #undef CANCEL
}


/**
    Optimizers for problems dealing with weighted 2 dimensional versions of TV norm regularization.
    One dimensional weighted solvers are used as building blocks within a proximal stacking framework.
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
void DR_proxDiff(size_t n, double* input, double* output, double* W, Workspace *ws);
void DR_columnsPass(size_t M, size_t N, double* input, double* output, double* W, Workspace **ws);
void DR_rowsPass(size_t M, size_t N, double* input, double* output, double* ref, double* W, Workspace **ws);

/* DR2L1W_TV
 *
 * This code puts in vector s the solution w for the problem
 * min_w gam(f_1(w)+f_2(w)) + 0.5*norm(w)^2 -u^w, where
 * f_1 and f_2 are cut functions (i.e., weighted TV terms) over cols and rows of 'y'
 * Thus, in other words, it solves a rewritten version of a weighted 2D-TV-L1!
 *
 * Inputs:
 * int M -- num of rows in input matrix (unary)
 * int N -- num of cols in input matrix
 * double* unary -- input matrix / image
 * double *W1    -- regularization weight along cols, as a matrix of size (M-1)xN
 * double *W2    -- regularization weight along rows, as a matrix of size Mx(N-1)
 * int nThreads  -- num of threads to use for parallel processing
 *
 * Outputs:
 *   s -- array containing primal solution (pre-alloced space assumed)
 */

int DR2L1W_TV(size_t M, size_t N, double*unary, double*W1, double*W2, double*s, int nThreads, int maxit, double* info)
{

  int i;
  double *ytr = NULL;
  double *t = NULL;
  double *tb = NULL;
  Workspace **ws = NULL;
  int maxDim;
  
  #define FREE \
    if(t) free(t); \
    if(tb) free(tb); \
    if(ws) freeWorkspaces(ws,nThreads);
    
  #define CANCEL(txt,info) \
      printf("DR2L1W_TV: %s\n",txt); \
      FREE \
      if(info) info[INFO_RC] = RC_ERROR;\
      return 0;

  //printMatrix(M, N, unary, "unary");

  if (nThreads < 1) nThreads = 1;
  omp_set_num_threads(nThreads);
  maxDim = (M > N) ? M : N;

  // Alloc memory for algorithm */
  t = (double*) malloc(sizeof(double)*M*N);
  tb = (double*)malloc(sizeof(double)*M*N);
  ws = newWorkspaces(maxDim,nThreads);

  if (!t || !tb || !ws)
    {CANCEL("out of memory", info)}
    
  /* Set number of iterations */
  if(maxit <= 0) maxit = MAX_ITERS_DR;

  // t = ones(size(unary)) * mean(unary(:));
  double sum=0;
  for (i=0; i < M*N; i++)
    sum += unary[i];
  sum = 2*sum / (M*N); 
  for (i=0; i < M*N; i++)
    t[i]=sum;

    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Starting Douglas-Rachford with size=[%d,%d], norms=[%lf,%lf], threads=%d\n",M,N,1,1,nThreads); fflush(DEBUG_FILE);
    #endif

  int iter = 0;
  /* MAIN LOOP */
  while ( iter < maxit ) {
    iter++;
    // reflect for B_vertical
    //s = 2*dualprojLines(t, u1, W1) - t;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Dual projection along columns\n"); fflush(DEBUG_FILE);
    #endif
    // Projection (prox) step
    DR_columnsPass(M, N, t, s, W1, ws);
    // Reflection
    for (i=0; i < M*N; i++) s[i] = 2*s[i] - t[i];
    
    // reflect for -B_horizontal 
    // t = 2*( -dualprojLines(-s', u2', W2')) - s';
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Dual projection along rows\n"); fflush(DEBUG_FILE);
    #endif
    // Projection (prox) step, taking into account displacemente from reference unary signal
    DR_rowsPass(M, N, s, tb, unary, W2, ws);
    // Reflection
    for (i=0; i < M*N; i++) tb[i] = -2*tb[i] - s[i];

    // Combiner step
    for (i=0; i < M*N; i++) t[i] = 0.5*(t[i]+tb[i]);
  }
  
  // DR is divergent, but with an additional projection we can recover valid solutions
  DR_columnsPass(M, N, t, s, W1, ws);
  DR_rowsPass(M, N, s, tb, unary, W2, ws);
  for (i = 0; i < M*N; i++) s[i] = - s[i] - tb[i];
  
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
    Performs the columns proximity step in the DR algorithm.

    @param M number of rows in input signal
    @param N number of columns in input signal
    @param input signal
    @param output signal
    @param W regularization weight along cols, as a matrix of size (M-1)xN
    @param ws array of Workspaces to use for the computation
*/
void DR_columnsPass(size_t M, size_t N, double* input, double* output, double* W, Workspace **ws) {
    #pragma omp parallel shared(M,N,input,output,W,ws) default(none) 
    {
        int i,j;
        // Get thread number
        int id = omp_get_thread_num();
        // Get corresponding workspace
        Workspace *wsi = ws[id];
        wsi->warm = 0;
        
        // Run 1-d solvers in parallel on each column of the input
        #pragma omp for
        for (j=0; j < N; j++) {
            resetWorkspace(wsi);
            // Array for weights
            double* wline = getDoubleWorkspace(wsi);
            // Prepare weights
            memcpy(wline, W+((M-1)*j), sizeof(double)*(M-1));
            // Prepare inputs
            memcpy(wsi->in, input+(M*j), sizeof(double)*M);
            // Compute prox difference for this column
            DR_proxDiff(M, wsi->in, wsi->out, wline, wsi);
            // Save output
            memcpy(output+(M*j), wsi->out, sizeof(double)*M);
        }
    }
}

/** 
    Performs the rows proximity step in the DR algorithm.

    @param M number of rows in input signal
    @param N number of columns in input signal
    @param input signal
    @param output signal
    @param ref reference signal (input to the global DR algorithm)
    @param W regularization weight along rows, as a matrix of size Mx(N-1)
    @param ws array of Workspaces to use for the computation
*/
void DR_rowsPass(size_t M, size_t N, double* input, double* output, double* ref, double* W, Workspace **ws) {
    #pragma omp parallel shared(M,N,input,ref,output,W,ws) default(none) 
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
            // Array for weights
            double* wline = getDoubleWorkspace(wsi);
            // Prepare weights
            int idx;
            for ( idx = j, i = 0 ; i < N-1 ; i++, idx+=M )
                wline[i] = W[idx];
            // Prepare inputs, considering displacement from reference signal
            for ( idx = j, i = 0 ; i < N ; i++, idx+=M )
                wsi->in[i] = ref[idx]-input[idx];
            // Compute prox difference for this row
            DR_proxDiff(N, wsi->in, wsi->out, wline, wsi);
            // Save output, recovering displacement from reference signal
            for ( idx = j, i = 0 ; i < N ; i++, idx+=M )
                output[idx] = wsi->out[i] - ref[idx];
        }
    }
}

/**
 Applies the weighted L1-TV proximity operator to the given 1D input and returns the difference between input and output of prox.
 
 @param n length of input signal
 @param input signal
 @param output signal
 @param W weights of the TV regularization
 @param ws Workspace to use for the computation
 */
void DR_proxDiff(size_t n, double* input, double* output, double* W, Workspace *ws) {
    int i;

    // Compute proximity
    tautString_TV1_Weighted(input, W, output, n);
    // Return differences between input and proximity output
    for (i=0; i < n; i++)
      output[i] = input[i] - output[i];   
}


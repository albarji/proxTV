#ifndef _TVOPT_CLEAN_H
#define _TVOPT_CLEAN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "general.h"
#include "LPopt.h"
#include "utils.h"
#include "condat_fast_tv.h"



/* Uncomment to print debug messages to a debug file */
//#define DEBUG
/* Choose here the name of the debug file */
#ifdef DEBUG
    static FILE* DEBUG_FILE = fopen("debug.tmp","w");
#else
    static FILE* DEBUG_FILE = NULL;
#endif
#define DEBUG_N 10 /* Maximum vector length to print in debug messages */

/* Uncomment to produce timing messages (in some of the algorithms) */
//#define TIMING

/* Includes for parallel computation */
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_set_num_threads(nThreads) /**/;
#endif

/* Indexes of info structure of solvers */
#define N_INFO 3
#define INFO_ITERS 0
#define INFO_GAP 1
#define INFO_RC 2

/* Comparison tolerance */
#define EPSILON 1e-10
#define IS_ZERO(x) (x < EPSILON & x > -EPSILON) 

/* Return Codes */
#define RC_OK 0 // Solution found at the specified error level
#define RC_ITERS 1 // Maximum number of iterations reaches, result might be suboptimal
#define RC_STUCK 2 // Algorithm stuck, impossible to improve the objetive function further, result might be suboptimal
#define RC_ERROR 3 // Fatal error during the algorithm, value of the solution is undefined.

#endif

/* TV-L1 solver */

/* Stopping tolerance */
#define STOP_PN 1e-6
/* Minimum decrease */
#define SIGMA 0.0500
/* Maximum number of iterations */
#define MAX_ITERS_PN 100

/* TV-L2 solver */

/* Stopping tolerances */
#define STOP_MS 1e-5 // Duality gap
#define STOP_MSSUB 1e-6 // Distance to boundary (prerequisite for termination)
/* Maximum number of iterations */
#define MAX_ITERS_MS 100

/* TV-Lp solver */

/* Stopping tolerance */
#define STOP_TVLP 1e-5 // Duality gap
/* Maximum number of iterations */
#define MAX_ITERS_TVLP 10000
#define MAX_ITERS_TVLPFW 1000000
#define MAX_ITERS_TVLPGPFW 1000000
/* Maximum number of iterations without improvement */
#define MAX_NOIMP_TVLP 10
/* Accuracy of internal prox LP solver for GP algorithm */
#define OBJGAP_LPPROX_TVLP 1e-15
/* Lambda approximation steps */
#define LAMBDA_STEPS_TVLP 1
/* Reduction factor for initial lambda approximation step */
#define LAMBDA_REDUCTION_TVLP 1e3
/* Stopping tolerance in stepsize magnitude for FW algorithm */
#define STOP_STEP_TVLP_FW 1e-15
/* Numer of FW cycles per GP step in hybrid GP+FW method */
#define FW_CYCLES_TVLP 10
/* Maximum number of iterations without significant improvement in hybrid GP+FW method */
#define MAX_NOIMP_TVLP_GPFW (10*FW_CYCLES_TVLP)
/* Minimum improvement to be considered significant */
#define MIN_IMP_TVLP 1e-10
/* Enforced GP cycles in GP+FW method when ill-conditioning is detected */
#define ENFORCED_CYCLES_TVLP_GPFW 10

/* Multidimensional TV solver */

/* Stopping tolerance */
#define STOP_PD 1e-6 // Mean absolute change in solution
/* Maximum number of iterations */
#define MAX_ITERS_PD 100
/* Maximum number of iterations for inner solvers */
#define MAX_ITERS_PD_INNER 0 //0

/* DR parameters */
#define STOP_DR 1e-6
#define MAX_ITERS_DR 100


/*** Function headers ***/

/* TV-L1 solvers */
int PN_TV1(double *y,double lambda,double *x,double *info,int n,double sigma,Workspace *ws);
int PN_TV1_Weighted(double* Y, double* W, double* X, double* info, int n, double sigma, Workspace* ws);
/* TV-L2 solvers */
int more_TV2(double *y,double lambda,double *x,double *info,int n);
int morePG_TV2(double *y,double lambda,double *x,double *info,int n,Workspace *ws);
int PG_TV2(double *y,double lambda,double *x,double *info,int n);
/* TV-Lp solvers */
int GP_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);
int OGP_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);
int FISTA_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);
int FW_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);
int GPFW_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);
/* General-dimension TV solvers */
int PD_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters);
int PD2_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters);
int DR_TV(double *y,double *lambdas,double *norms,double *dims,double *x,double *info,int *ns,int nds,int npen,int ncores,int maxIters);
int proxDykstraTV2DWeighted(double* y, double* W, double* x,int* ns, int nds, double* info,int ncores,int maxIters);

#endif

/**
    General definitions
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#ifndef _GENERAL_H
#define _GENERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/* Mex and LAPACK includes */
#ifdef NOMATLAB
#undef lapack_int
#define lapack_int              long
#include <lapacke.h>
inline double mxGetInf() { return INFINITY; }
#else
#include "mex.h"
#include "lapack.h"
#include "matrix.h"
#define lapack_int  ptrdiff_t
#endif

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
#define IS_POSITIVE(x) (x > EPSILON)
#define IS_NEGATIVE(x) (x < -EPSILON)

/* Return Codes */
#define RC_OK 0 // Solution found at the specified error level
#define RC_ITERS 1 // Maximum number of iterations reaches, result might be suboptimal
#define RC_STUCK 2 // Algorithm stuck, impossible to improve the objetive function further, result might be suboptimal
#define RC_ERROR 3 // Fatal error during the algorithm, value of the solution is undefined.

#endif

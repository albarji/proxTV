/* File: proxtv.cc -*- c++ -*- */
/* (c) Copyright 2013  Alvaro Barbero and Suvrit Sra */
/* C++ library header */

#include "proxtv.h"

/* Solves: 0.5*norm(x-y)^2 + lam*TVL1(x)
 * 
 * INPUT
 * y --   vector of doubles
 * n --   length of y
 * lam -- regularization parameter
 *
 * OUTPUT
 * x -- solution vector
 * info -- status info about solver / solution
 *
 * ASSUMPTION
 * x -- preallocated array of length 'n'
 * 
 */
int tv1d (double* y, double* x, size_t n, double lam, double p, double* info)
{
  if (p==1) {
    TV1D_denoise(y, x, n, lam);
    info=0;
    return 0;
  } else if (p==2) { 
    return morePG_TV2(y, lam, x, info, n, NULL);
  } else {
    return GPFW_TVp(y, lam, x, info, n, p, NULL);
  }
}


/* Solves weighted-1D-TV
   0.5*norm(x-y)^2 + sum_i lam(i)*|x(i+1)-x(i)|
 * 
 * INPUT
 * y --   vector of doubles
 * n --   length of y
 * lam --   vector of weights parameter
 *
 * OUTPUT
 * x -- solution vector
 * info -- status info about solver / solution
 *
 * ASSUMPTION
 * x -- preallocated array of length 'n'
 * 
 */
int tv1dw(double* y, double* x, size_t n, double* lam, double* info)
{
  return PN_TV1_Weighted(y, lam, x, info, n, SIGMA, NULL);
}

/* Solves 2D-TV
 * 0.5*norm(x-y)^2 + sum_i lam*|x(i+1,j)-x(i,j)| + sum_j lam*|x(i,j+1)-x(i,j)|
 * 
 * INPUT
 * y --   vector of doubles
 * n --   length of y
 * lam -- regularization parameter
 * nthr --num of threads to use
 * maxit--max num of iterations of underlying solver
 *
 * OUTPUT
 * x -- solution vector
 * info -- status info about solver / solution
 *
 * ASSUMPTION
 * x -- preallocated array of length 'n'
 * 
 */
int tv2d (double* y, double* x, size_t m, size_t n, double* lam, double pr,
          double pc, int nthr, int maxit, double* info)
{
  return 0;
}

int tv2dw(double* y, size_t m, size_t n, double* lam, double pr,
          double pc, int nthr, int maxit, double* info)
{
  return 0;
}

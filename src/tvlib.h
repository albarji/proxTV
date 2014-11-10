/* File: tvlib.h -*- c++ -*- */
/* (c) Copyright 2013  Alvaro Barbero and Suvrit Sra */
/* C++ library header */

#ifndef _TVLIB_H_
#define _TVLIB_H_

#ifndef NOMATLAB
#define NOMATLAB
#endif 

#include "TVopt.h"

#define _TV_DEF_NTH 1
#define _TV_DEF_MIT 50

int tv1d (double* y, double*x, size_t n, double lam, double p=1, int nthr = _TV_DEF_NTH, int maxit=_TV_DEF_MIT, double* info=0);
int tv1dw(double* y, double*x, size_t n, double* lam, double p=1, int nthr = _TV_DEF_NTH, int maxit=_TV_DEF_MIT, double* info=0);
int tv2d (double* y, double*x, size_t m, size_t n, double* lam, double pr=1, double pc=1, int nthr = _TV_DEF_NTH, int maxit=_TV_DEF_MIT, double* info=0);
int tv2dw(double* y, double*x, size_t m, size_t n, double* lam, double pr=1, double pc=1, int nthr = _TV_DEF_NTH, int maxit=_TV_DEF_MIT, double* info=0);

#endif // _TVLIB_H_

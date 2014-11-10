/* File : proxtv.i */
%module proxtv_internal
%include "carrays.i"
%array_class(double, doubleArray);
%{
#define SWIG_FILE_WITH_INIT
#include "proxtv.h"
extern int tv1d (double* y, double* x, size_t n, double lam, double p, double* info);
extern int tv1dw(double* y, double* x, size_t n, double* lam, double* info);
%}

extern int tv1d (double* y, double* x, size_t n, double lam, double p, double* info);
extern int tv1dw(double* y, double* x, size_t n, double* lam, double* info);

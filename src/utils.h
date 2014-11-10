/**
    General utility functions.
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/*** Structures ***/

/* Workspace for warm restarts and memory management */
typedef struct {
    /* Size of memory vectors */
    int n;
    /* Generic memory which can be used by 1D algorithms */
    double **d;
    int maxd, nd;
    int **i;
    int maxi, ni;
    /* Memory for inputs and outputs */
    double *in,*out;
    /* Warm restart variables */
    short warm;
    double *warmDual;
    double warmLambda;
} Workspace;

/* Workspace defines */
#define WS_MAX_MEMORIES 100

/*** Module functions ***/

/* Math functions */
short sign(double s);
double min(double x, double y);
double max(double x, double y);
void radialProjection(double *x, int n, double norm, double lambda);

/* Memory allocation functions */
Workspace* newWorkspace(int n);
void resetWorkspace(Workspace *ws);
double* getDoubleWorkspace(Workspace *ws);
int* getIntWorkspace(Workspace *ws);
void freeWorkspace(Workspace *ws);
Workspace** newWorkspaces(int n,int p);
void freeWorkspaces(Workspace **wa,int p);

/* Other utility functions */
int compareDoublesReversed(const void *v1, const void *v2);

#endif

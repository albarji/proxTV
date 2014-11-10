/**
    Optimizers for problems dealing with Lp norm regularization or Lp ball constraints.
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#ifndef _LPOPT_H
#define _LPOPT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "general.h"
#include "utils.h"

/*** Definitions ***/

/* Lp solver */

/* Stopping tolerance */
#define STOP_PNLP 0 //1e-12 //Relative change in objective function
#define STOP_GAP_PNLP 1e-5 //Absolute dual gap for stopping
/* Maximum number of iterations */
#define MAX_ITERS_PNLP 1000 //100
/* Threshold for Armijo sufficient descent condition in stepsize selection */
#define SIGMA_PNLP 0.05
/* Tolerance for considering an x component as 0 */
#define EPSILON_PNLP 1e-15
/* Minimum acceptable gradient * direction value for Hessian directions */
#define MIN_GRD_HESSIAN_PNLP 1e-15
/* Approximation limits for Lp norm */
#define LPPROJ_PSMALL 1.002
#define LPPROJ_PLARGE 100
/* Minimum stepsize */
#define MIN_STEP_PNLP 1e-10

/*** Module functions ***/

/* Lp norm value */
double LPnorm(double *x, int n, double p);

/* Lp proximity solver */
int PN_LP1(double *y,double lambda,double *x,double *info,int n);
int PN_LP2(double *y,double lambda,double *x,double *info,int n);
int PN_LPinf(double *y,double lambda,double *x,double *info,int n,Workspace *ws);
int PN_LPp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws,int positive,double objGap);
int PN_LPp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws,int positive);
/* Lp-ball projections */
int LP1_project(double *y,double lambda,double *x,int n,Workspace *ws);
int LPp_project(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws);

/* Linear LP constrained solver */
void solveLinearLP(double *z, int n, double p, double lambda, double *s);

#endif

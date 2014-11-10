/**
    Optimizers for problems dealing with Lp norm regularization or Lp ball constraints.
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "LPopt.h"

/*** Internal functions headers ***/
double PN_LPpGap(double *x, double *y, double *diff, int n, double q, double lambda, double norm);

/*** Internal definitions ***/

/* Kind of PNLP update */
#define PNLP_HESSIAN 0 // Hessian update
#define PNLP_MNSG 1 // Minimum Norm SubGradient update
#define PNLP_GRAD 2 // Gradient update
#define PNLP_BB 3 // Barzilai-Borwein update

/*** Module functions ***/

/**
    Computes the value of the Lp norm for a given vector.
    If the value of p is too similar to 1, the norm is approximated as an L1 norm.
    If the value of p is too large, the norm is approximated as an Linf norm.
    For the general case, instead of computing the norm using the standard formula, a normalized version is used.
    This ensures a more stable computation. The formula used is:
    
        norm(x,p) = norm(x,inf) * (sum_i |x(i) / norm(x,inf)|^p)^(1/p)
        
    Inputs:
        - double *x: vector for which to compute the norm.
        - double p: Lp norm value.
        
    Output: Lp norm value of x.
*/
double LPnorm(double *x, int n, double p) {
    double norm = 0;
    
    /* Small norm limit case */
    if ( p <= LPPROJ_PSMALL ) {
        for ( int i = 0 ; i < n ; i++ )
            norm += fabs(x[i]);
        return norm;    
    }
    
    /* Compute normalization factor c = norm(x,inf) = max_i |x(i)| */
    double aux, c=0;
    for ( int i = 0 ; i < n ; i++ ) {
        aux = fabs(x[i]);
        if ( aux > c )
            c = aux;
    }
    
    /* If c == 0, then the norm is 0, regardless of p */
    if ( c == 0 )
        return 0;
        
    /* Large norm limit case: just return c */
    if ( p >= LPPROJ_PLARGE ) {
        return c;
    }
   
    /* General case: normalized norm computation */
    for ( int i = 0 ; i < n ; i++ )
        norm += pow( fabs(x[i]/c) , p );
    norm = c*pow( norm , 1./p );
    
    return norm;
}

/*  PN_LP1

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x||_1 ,
        
    for the L1 norm. This is done in closed form using the soft-thresholding operator
    
        x = sgn(y) .* max {|y| − lambda, 0}
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
*/
int PN_LP1(double *y,double lambda,double *x,double *info,int n){
    /* Apply soft-thresholding iteratively */
    for (int i = 0 ; i < n ; i++ )
        x[i] = sign(y[i]) * max(fabs(y[i]) - lambda, 0);
    
    /* Fill info fields */
    if (info) {
        info[INFO_ITERS] = 0;
        info[INFO_GAP] = 0;
        info[INFO_RC] = RC_OK;
    }
    
    return 1;
}

/*  PN_LP2

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x||_2 ,
        
    for the L2 norm. This is done in closed form using the operator
    
        x = y * max {||y||_2 - lambda , 0} / ||y||_2
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
*/
int PN_LP2(double *y,double lambda,double *x,double *info,int n){
    /* Compute norm of y */
    double norm = LPnorm(y, n, 2);
    /* Special case where norm = 0 : solution is an all zeros vector */
    if ( norm == 0 ) {
        memset(x, 0, sizeof(double)*n);
    /* General case */
    } else {
        /* Correction factor */
        double corr = max( norm - lambda , 0 );

        /* Apply close-form solution operator for each entry */
        for (int i = 0 ; i < n ; i++ )
            x[i] = y[i] * corr / norm;
    }
    
    /* Fill info fields */
    if (info) {
        info[INFO_ITERS] = 0;
        info[INFO_GAP] = 0;
        info[INFO_RC] = RC_OK;
    }
    
    return 1;
}

/*  PN_LPinf

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x||_inf ,
        
    for the Linf norm. This is done by means of solving the easier dual projection problem
    
        min_z 0.5 ||z-y||^2  s.t. ||z||_1 <= lambda
        
    The primal solution is then recovered using Moreu's decomposition
    
        x = y - z
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
*/
int PN_LPinf(double *y,double lambda,double *x,double *info,int n,Workspace *ws){
    /* Invoke dual L1 projection problem */
    LP1_project(y, lambda, x, n, ws);
    
    /* Recover primal solution using Moreau's decomposition */
    for(int i = 0 ;i < n ; i++)
        x[i] = y[i] - x[i];
    
    /* Fill info fields */
    if (info) {
        info[INFO_ITERS] = 0;
        info[INFO_GAP] = 0;
        info[INFO_RC] = RC_OK;
    }
    
    return 1;
}

/*  PN_LPp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x||_p ,
        
    for the Lp norm. To do so a Projected Newton algorithm is used.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the Lp norm.
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
        - positive: 1 if all inputs y >= 0, 0 else.
        - objGap: desired quality of the solution in terms of duality gap.
*/
int PN_LPp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws,int positive,double objGap){
    double *g=NULL,*d=NULL,*xnorm=NULL,*auxv=NULL;
    double stop,stop2,q,qInv,nx,f,fupdate,aux,c,den,xp1vGrad,gRd,delta,prevDelta,improve,rhs,grad0,gap,epsilon;
    int *inactive=NULL,*signs=NULL;
    int i,j,iters,recomp,found,nI;
    short updateKind,stuck;
    
    #define CHECK_INACTIVE(x,g,inactive,nI,i) \
        for(i=nI=0 ; i<n ; i++) \
            if( x[i] > epsilon || (x[i] <= epsilon && g[i] < -EPSILON) )  \
                inactive[nI++] = i;
        
    #define FREE \
        if(!ws){ \
            if(g) free(g); \
            if(d) free(d); \
            if(xnorm) free(xnorm); \
            if(auxv) free(auxv); \
            if(inactive) free(inactive); \
            if(signs) free(signs); \
        }
        
    #define CANCEL(txt,info) \
        printf("PN_LPp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Reduction of stepsize if no improvement detected */
    #define NO_IMPROVE_CUT 0.1
    
    /* Compute dual norm q */
    q = 1/(1-1/p);
    qInv = -1/q;
    
    /* Special case where the solution is the trivial x = 0 */
    /* This is bound to happen if ||y||_q <= lambda */
    if ( LPnorm(y, n, q) <= lambda ) {
        /* Set solution to all zeros */
        memset(x, 0, sizeof(double)*n);
        /* Fill return info */
        if (info) {
            info[INFO_RC] = RC_OK;
            info[INFO_GAP] = 0;
            info[INFO_ITERS] = 0;
        }
        FREE
        return 1;
    }

    /* Alloc memory if no workspace available */
    if(!ws){
        g = (double*)malloc(sizeof(double)*n);
        d = (double*)malloc(sizeof(double)*n);
        xnorm = (double*)malloc(sizeof(double)*n);
        auxv = (double*)malloc(sizeof(double)*n);
        inactive = (int*)malloc(sizeof(int)*n);
        if(!positive){
            signs = (int*)malloc(sizeof(int)*n);
            if(!signs){CANCEL("out of memory",info)}
        }
    }
    /* If a workspace is available, assign pointers */
    else{
        g = getDoubleWorkspace(ws);
        d = getDoubleWorkspace(ws);
        xnorm = getDoubleWorkspace(ws);
        auxv = getDoubleWorkspace(ws);
        inactive = getIntWorkspace(ws);
        if(!positive) {
            signs = getIntWorkspace(ws);
            if(!signs){CANCEL("out of memory",info)}
        }
    }
    if(!g || !d || !xnorm || !auxv || !inactive)
        {CANCEL("out of memory",info)}
        
    /* As initial points, use an approximation for the closed-form solution */

    /* Find close-form solution for TV-L2 */
    PN_LP2(y , lambda , x, NULL ,n);
    /* Measure gap of this solution in terms of TV-Lp norm */
    for ( i = 0 ; i < n ; i++ )
        auxv[i] = x[i] - y[i];
    nx = LPnorm( x, n, p );
    stop = PN_LPpGap(x, y, auxv, n, q, lambda, nx);
    
    /* Find close-form solution for TV-L1 or TV-Linf, whichever is closer in norm value */
    if ( p < 2 )
        PN_LP1(y , lambda , xnorm, NULL ,n);
    else
        PN_LPinf(y , lambda , xnorm, NULL ,n, ws);
    /* Measure gap of this solution in terms of TV-Lp norm */
    for ( i = 0 ; i < n ; i++ )
        auxv[i] = xnorm[i] - y[i];
    nx = LPnorm( xnorm, n, p );
    stop2 = PN_LPpGap(xnorm, y, auxv, n, q, lambda, nx);
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Starting candidates: p2 = %lf, p%s = %lf\n", stop, p<2 ? "1" : "inf", stop2);
    #endif
    
    /* Use as starting point the closed-form solution with smaller dual gap */
    if ( stop2 < stop ) {
        /* Copy starting point */
        memcpy(x, xnorm, sizeof(double)*n);
        stop = stop2;
    }
    /* If dual gap is good enough, we are finished */
    if ( stop < objGap ) {
        FREE
        if (info) {
            info[INFO_GAP] = stop;
            info[INFO_ITERS] = 0;
        }
        return 1;
    }
    
    /* Take out the sign from the input, remembering original signs */
    /* Remove also from the starting point */
    if(!positive){
        for(i=0;i<n;i++){
            signs[i] = sign(y[i]);
            x[i] = fabs(x[i]);
            y[i] = fabs(y[i]);
        }
    }
    
    /* Tolerance for considering an x component as 0 */
    epsilon = EPSILON_PNLP;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"p=%lf, q=%lf, qInv=%lf\n",p,q,qInv);
        fprintf(DEBUG_FILE,"lambda=%lf\n",lambda);
        fprintf(DEBUG_FILE,"epsilon=%lg\n",epsilon);
        fprintf(DEBUG_FILE,"y=[ ",iters);
        for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lf ",y[i]);
        fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"Start, x=[ ",iters);
        for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lf ",x[i]);
        fprintf(DEBUG_FILE,"]\n");
    #endif
    
    /* Initial check for zero values */
    for ( i = 0 ; i < n ; i++ )
        if ( x[i] < epsilon )
            x[i] = epsilon;
    
    /* Initial value of the point norm */
    nx = LPnorm(x,n,p);
    /* Compute differences x - y */
    for ( i = 0 ; i < n ; i++ )
        auxv[i] = x[i]-y[i];  /* auxv storing differences x-y*/
        
    /* Projected Newton loop */
    stop = gap = DBL_MAX; iters = 0;
    for(iters=0 ; stop > STOP_PNLP && iters < MAX_ITERS_PNLP && gap > objGap ; iters++){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, x=[ ",iters);
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",x[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Compute current obj. function value */

        /* Current obj. function value */
        f = 0;
        for(i=0;i<n;i++)
            f += auxv[i]*auxv[i];
        f = 0.5*f + lambda*nx;

        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, fval=%lf\n",iters,f);
        #endif
        
        /* Gradient at the current x */
        for(i=0;i<n;i++){
            // Compute normalized current vector, powered to p-1
            xnorm[i] = pow(x[i] / nx , p-1);
            // Compute gradient
            g[i] = auxv[i] + lambda * xnorm[i];
        }

        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, xnorm=[ ",iters);
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",xnorm[i]);
            fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"Iter %d, g=[ ",iters);
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",g[i]);
            fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"||g||_inf = %lg, ||g||_2=%lf\n",LPnorm(g,n,DBL_MAX),LPnorm(g,n,2));
        #endif
        
        #ifdef DEBUG
        {
            double drift = 0;
            for(i=0;i<n;i++)
                drift += fabs(g[i] - (x[i] - y[i] + lambda * pow(x[i] / nx , p-1)));
            fprintf(DEBUG_FILE,"Iter %d, gDrift=%lf\n",iters,drift);
        }
        #endif
        
        /* Find inactive constraints (variables to update) */
        CHECK_INACTIVE(x,g,inactive,nI,i)        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, nI=%d\n",iters,nI);
        #endif
        /* Special case where no inactive constraints exist: the optimum has been found */
        if(!nI) break;
        
        /* If the norm of x is zero, the Hessian is undefined */
        for ( i = 0 ; i < nI ; i++ ) {
            j = inactive[i];
            if ( x[j] > epsilon )
                break;
        }
        /* In that case, enforce minimum norm subgradient update */
        if ( i == nI ) {
            updateKind = PNLP_MNSG;
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"||x|| = 0, enforcing minumum norm subgradient step\n");
            #endif
        }
        else {
            updateKind = PNLP_HESSIAN;
        }
            
        /* Calculate updating direction using the chosen strategy */
        switch (updateKind) {
        
            /* Use the gradient as the (reversed) updating direction */
            case PNLP_GRAD: 
                /* Copy over the gradient for inactive variables */
                for(i=0;i<nI;i++){
                    j = inactive[i];
                    d[j] = g[j];
                }
                break;
                
            /* Use minimum norm subgradient updating direction */
            case PNLP_MNSG:
                /* Compute minimum norm subgradient at x=0 */
                /* Since the gradient at x=0 is not to be trusted, the inactive constraint detection is ignored */
                nI = n;
                for(i=0;i<n;i++){
                    d[i] = - pow(y[i] / lambda, 1 / (p-1));
                    inactive[i] = i;
                }
                break;
            
            /* Hessian (Newton) updating direction */
            case PNLP_HESSIAN:
                /* Auxiliar terms for updating direction */
                c = lambda * (1-p) * 1./nx;
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"Iter %d, c=%lg\n",iters,c);
                #endif
            
                /* Precompute auxv = 1./(1-c * xnorm./x) .* xnorm vector */
                /* At the same time, compute the 1./(1-c xnorm./x) .* g vector, which is the diagonal part of the direction
                   and the inner product  1./(1-c * xnorm./x) .* xnorm · g */
                xp1vGrad = 0;
                for ( i = 0 ; i < nI ; i++ ) {
                    j = inactive[i];
                    /* Common component */
                    auxv[j] = 1. / (1. - c * pow(x[j] / nx , p-2) );
                    /* 1./(1+c xnorm./x) .* g vector */
                    d[j] = auxv[j] * g[j];
                    /* 1./(1+c * xnorm./x) .* xnorm */
                    auxv[j] *= xnorm[j];
                    /* 1./(1+c * xnorm./x) .* xnorm · g */
                    xp1vGrad += auxv[j] * g[j];
                }
                
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"Iter %d, dDiag=[ ",iters);
                    for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",d[i]);
                    fprintf(DEBUG_FILE,"]\n");
                    fprintf(DEBUG_FILE,"Iter %d, v·xnorm=[ ",iters);
                    for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",auxv[i]);
                    fprintf(DEBUG_FILE,"]\n");
                    fprintf(DEBUG_FILE,"xp1vGrad=%lg\n",xp1vGrad);
                #endif
                
                /* Precompute denominator of 1-rank part of direction */
                den = 1./c;
                for(i=0;i<nI;i++){
                    j = inactive[i];
                    den += xnorm[j] * auxv[j];
                }
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"Iter %d, den=%lg\n",iters,den);
                #endif
                
                /* 1-rank part of direction */
                den = xp1vGrad / den;
                for(i=0;i<nI;i++){
                    j = inactive[i];
                    /* For auxv = 0 components, the Hessian is unstable. Instead, use the gradient */
                    if ( ! auxv[j] )
                        d[j] = g[j];
                    /* Else, use full Hessian direction */
                    else
                        d[j] -= den * auxv[j];
                }
            break;
        }
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, d=[ ",iters);
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",d[i]);
            fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"||d||_inf = %lg, ||d||_2=%lf\n",LPnorm(d,n,DBL_MAX),LPnorm(d,n,2));
        #endif
        
        /*** Stepsize computation ***/
        
        /* Gradient times direction */
        gRd = 0;
        /* In the MNSG update the gradient is not to be trusted, also, we want to get out of the current point as soon as possible,
           so we leave gRd = 0 to accept any improving step, no matter how small. The same happens for the BB iteration */
        /* In other updates the gradient is safe, so we resort to the usual calculation */
        if ( updateKind != PNLP_MNSG && updateKind != PNLP_BB ) {
            for(i=0;i<nI;i++){
                j = inactive[i];
                gRd += g[j] * d[j];
            }
        }
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"gRd=%lg\n", gRd);
        #endif
        
        /* Copy over the full current point, in order to have in auxv those components that wont be updated */
        for(i=0;i<n;i++)
            auxv[i] = x[i]; /* auxv will hold the tentative new point x + delta · d*/

        /* Iterate Armijo rule to find acceptable stepsize */
        recomp = 0; delta = 1; found = 0;
        while(!found){
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"Iter %d, delta=%lg\n",iters,delta);
            #endif
            
            /* Compute projected point after update for the current stepsize */
            for(i=0;i<nI;i++){
                j = inactive[i];
                auxv[j] = x[j] - delta * d[j];
                if(auxv[j] < epsilon) auxv[j] = epsilon;
            }
            
            /* Compute Lp norm after the update */
            nx = LPnorm(auxv,n,p);

            #ifdef DEBUG
                fprintf(DEBUG_FILE,"Iter %d, xP=[ ",iters);
                for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",auxv[i]);
                fprintf(DEBUG_FILE,"]\n");
                fprintf(DEBUG_FILE,"||xP||_%lg = %lg\n",p,nx); fflush(DEBUG_FILE);
            #endif
            
            /* If delta too small, the algorithm is unable to obtain further improvement */
            if(delta < MIN_STEP_PNLP) {
                /* Stop optimization and consider solution found */
                stop = 0;
                break;
            }

            /* Compute obj. value after the update */
            fupdate = 0;
            for(i=0;i<n;i++){
                aux = auxv[i]-y[i];
                fupdate += aux*aux;
            }
            fupdate = 0.5*fupdate + lambda*nx;
            
            /* Compute improvement */
            improve = f - fupdate;
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"f=%lg, fupdate=%lg, improve=%lg\n", f, fupdate, improve);
            #endif

            /* If negative improvement, force delta down and repeat */
            if(improve <= 0){
                delta *= NO_IMPROVE_CUT;
                continue;
            }
            
            /* Compute right hand side of Armijo rule */
            rhs = SIGMA_PNLP * delta * gRd;
            
            /* Check if the rule is met */

            #ifdef DEBUG
                fprintf(DEBUG_FILE,"rhs=%lg\n", rhs);
            #endif
            if(improve >= rhs) found = 1;
            else{
                /* If at start of step recomputation, compute gradient w.r.t delta at delta=0 */
                /* grad0 = -(xred - yred + lam * xnorm)'*d */
                if(!recomp){
                    grad0 = 0;
                    for ( i = 0 ; i < nI ; i++ ) {
                        j = inactive[i];
                        //grad0 -= (x[j] * y[j] + lambda * xnorm[j]) * d[j];
                        grad0 -= (x[j] - y[j] + lambda * xnorm[j]) * d[j];
                    }
                    recomp = 1;
                }
                /* Use quadratic interpolation to determine next stepsize */
                prevDelta = delta;
                aux = grad0 * delta;
                delta = - (aux*delta) / (2 * (-improve - aux));
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"Interpolated delta=%lg\n",delta);
                #endif
                /* If delta too similar to or larger than previous stepsize, or negative, halve it */
                if(prevDelta - delta < EPSILON || delta < 0) delta = prevDelta/2;
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"Corrected delta=%lg\n",delta);
                #endif
            }
        }
        
        /* Update point */
        for(i=0;i<nI;i++){
            j = inactive[i];
            x[j] = auxv[j];
        }
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, x=[ ",iters);
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",x[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Compute differences x - y of new point */
        /* This can be done efficiently by realizing that auxv already contains x */
        for ( i = 0 ; i < n ; i++ )
            auxv[i] -= y[i];
        
        /* Compute stopping criterion: relative change in objetive function */
        /* This is averaged with the stop value from the previous iteration,
           to avoid stopping at incidental low improvement iters */
        if ( stop < DBL_MAX )
            stop = 0.5 * improve / fupdate + 0.5 * stop;
        else
            stop = improve / fupdate;
        
        /* Compute dual gap */
        gap = PN_LPpGap(x, y, auxv, n, q, lambda, nx);
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, stop=%lg, gap=%lg\n",iters,stop,gap);
        #endif
    }
    
    /* Set almost zero entries to exact zero */
    for ( i = 0 ; i < n ; i++ )
        if (x[i]  <= epsilon)
            x[i] = 0;
    
    /* Place the sign back into the solution and the original data*/
    /* (only if input was not all positive) */
    if(!positive){
        for(i=0;i<n;i++){
            if(signs[i] == -1){ 
                x[i] = -x[i];
                y[i] = -y[i];
            }
        }
    }
    
    /* Check termination condition */
    if ( iters >= MAX_ITERS_PNLP ) {
        if (info) info[INFO_RC] = RC_ITERS;
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"WARNING: top limit of iterations reached\n"); fflush(DEBUG_FILE);
        #endif
    } else if ( gap > objGap ) {
        if (info) info[INFO_RC] = INFO_GAP;
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"WARNING: algorithm stuck, suboptimal solution found\n"); fflush(DEBUG_FILE);
        #endif   
    } else {
        if (info) info[INFO_RC] = RC_OK;
    }
    
    /* Free resources and return */
    FREE
    if (info) {
        info[INFO_GAP] = gap;
        info[INFO_ITERS] = iters;
    }
    return 1;
    
    #undef CHECK_INACTIVE
    #undef FREE
    #undef CANCEL
    #undef NO_IMPROVE_CUT
}

/*  PN_LPp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x||_p ,
        
    for the Lp norm. To do so a Projected Newton algorithm is used.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the Lp norm.
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
        - positive: 1 if all inputs y >= 0, 0 else.
        
    The proximity problem is solved to a default level of accuracy, as given by STOP_GAP_PNLP.
*/
int PN_LPp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws,int positive) {
    return PN_LPp(y, lambda, x, info, n, p, ws, positive, STOP_GAP_PNLP);
}

/** PN_LPpGap

    Computes the dual gap of the current PN-Lp prox solution.
    
    @argument x current PN-Lp prox solution
    @argument y reference point for PN-Lp problem
    @argument diff differences vector x - y
    @argument n length of x and y
    @argument q dual norm Lq
    @argument lambda multiplier of the Lp norm
    @argument norm precomputed Lp norm of x
**/
double PN_LPpGap(double *x, double *y, double *diff, int n, double q, double lambda, double norm) {
    /* Compute dual norm */
    double dualnorm = LPnorm(diff,n,q);
    
    /* Compute rescale factor to obtain feasible primal */
    double c;
    if ( dualnorm <= lambda )
        c = 1;
    else
        c = lambda / dualnorm;
        
    /* Initialize gap */
    double gap = lambda * norm;
        
    /* Add (1+c^2) 1/2 ||diff||^2 term */
    /* Also add c * y' diff term */
    double coef = (1.+c*c) * 0.5;
    for ( int i = 0 ; i < n ; i++ ) {
        gap += coef * diff[i] * diff[i];
        gap += c * y[i] * diff[i];
    }
    
    /*
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"(PN_LPpGap) norm=%lg, ||x||=%lg\n",norm,LPnorm(x,n,q/(q-1.0)));
        fprintf(DEBUG_FILE,"(PN_LPpGap) x= [ "); for(int i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",x[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"(PN_LPpGap) y= [ "); for(int i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",y[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"(PN_LPpGap) (x-y)= [ "); for(int i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%lg ",diff[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"(PN_LPpGap) c=%lg, coef=%lg\n",c,coef);
        fprintf(DEBUG_FILE,"(PN_LPpGap) gap=%lg\n",gap);
        fflush(DEBUG_FILE);
    #endif
    */
    
    return fabs(gap);
}

/*  LP1_project

    Given a point y, computes its euclidean projection onto the ||x||_1 = lambda ball
    
        min_x ||x-y||^2   s.t.  ||x||_1 <= lambda .
        
    This is solved efficiently using John Duchi's method. The code in this function is
    a C port of the original Matlab code.
    
    Inputs:
        - y: reference point.
        - lambda: ball radius.
        - x: array in which to store the projection.
        - n: length of array y (and x).
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
        
    Output:
        int with 1 value if solution found, 0 if error.
*/
int LP1_project(double *y,double lambda,double *x,int n,Workspace *ws){
    double *u=NULL,*sv=NULL;
    double theta;
    int i,rho;
    
    #define FREE \
        if(!ws){ \
            if(u) free(u); \
            if(sv) free(sv); \
        }
        
    #define CANCEL(txt) \
        printf("LP1_project: %s\n",txt); \
        return 0;
    
    // Lambda safety check
    if(lambda < 0) lambda = 0;
    
    // Take memory from workspace if available
    if(ws){
        u = getDoubleWorkspace(ws);
        sv = getDoubleWorkspace(ws);
    }
    else{
        u = (double*)malloc(sizeof(double)*n);
        sv = (double*)malloc(sizeof(double)*n);
    }    
    if(!u || !sv)
        {CANCEL("insufficient memory")}

    // Copy absolute of input values
    for(i=0;i<n;i++)
        u[i] = fabs(y[i]);
    // Sort them in descending order
    qsort((void*)u, n, sizeof(double),compareDoublesReversed);

    // Compute cumulative sum of sorted values
    sv[0] = u[0];
    for(i=1;i<n;i++)
        sv[i] = sv[i-1] + u[i];

    // Find thresholding index
    for(i=n-1,rho=-1 ; rho==-1 && i>=0 ; i--){
        if(u[i] > (sv[i] - lambda) / (i+1))
            rho = i;
    }

    // Compute thresholding value
    theta = (sv[rho] - lambda) / (rho+1);
    if(theta < 0) theta = 0;
    
    // Get projection
    for(i=0;i<n;i++){
        x[i] = fabs(y[i]) - theta;
        if(x[i] < 0) x[i] = 0;
        if(y[i] < 0) x[i] = -x[i];
    }
    
    // Free resources and return
    FREE
    return 1;
    
    #undef FREE
    #undef CANCEL
}

/*  LPp_project

    Given a point y, computes its euclidean projection onto the ||x||_p = lambda ball
    
        min_x ||x-y||^2   s.t.  ||x||_p <= lambda .
        
    This problem is dual to Lp proximity, and this duality is used to solve it through
    the Lp-prox solver.
    
    Inputs:
        - y: reference point.
        - lambda: ball radius.
        - x: array in which to store the projection.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the Lp ball.
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
*/
int LPp_project(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws){
    double q,norm,scale;
    int *s=NULL;
    int i;
    
    #define FREE \
        if(!ws) \
            if(s) free(s);
        
    #define CANCEL(txt,info) \
        printf("LPp_project: %s\n",txt); \
        *info = (double)RC_ERROR; \
        return 0;
        
    /* If p too small or too large, approximate by p=1 or p=inf */
    if(p <= LPPROJ_PSMALL) p = 1;
    else if(p >= LPPROJ_PLARGE) p = mxGetInf();
        
    /* Special case p == inf: solved by orthogonal projection in each component */
    if(p == mxGetInf()){
        /* Orthogonal projection */
        for(i=0;i<n;i++){
            if(y[i] > lambda) x[i] = lambda;
            else if(y[i] < -lambda) x[i] = -lambda;
            else x[i] = y[i];
        }
        return 1;
    }
    /* Special case p == 1: solved by using Duchi's method */
    else if(p == 1){
        if(!LP1_project(y,lambda,x,n,ws))
            {CANCEL("error in internal LP1 projection",info)}
        return 1;
    }
    /* In what follows: general p case */

    /* Compute dual norm */
    q = 1./(1.-1./p);
    
    /* Alloc memory for sign information (or take from workspace) */
    if(ws)
        s = getIntWorkspace(ws);
    else
        s = (int*)malloc(sizeof(int)*n);
    if(!s)
        {CANCEL("insufficient memory",info)}
    
    /* Remove sign from inputs */
    for(i=0;i<n;i++){
        s[i] = sign(y[i]);
        y[i] = fabs(y[i]);
    }
    
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(LPp_project), abs(y)=[ ");
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",y[i]);
            fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"(LPp_project), sign(y)=[ ");
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%d ",s[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif

    /* Invoke Lp prox solver on dual norm */
    if(!PN_LPp(y,lambda,x,info,n,q,ws,1))
        {CANCEL("error in internal Lp prox solver",info)}
    
    /* Apply Moreau's decomposition to recover primal problem solution */
    /* primal(y)* = y - dual(y)*     */
    for(i=0;i<n;i++)
        x[i] = y[i] - x[i];
        
    /* Place sign back */
    for(i=0;i<n;i++){
        if(s[i] < 0){
            y[i] = -y[i];
            x[i] = -x[i];
        }
    }
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(LPp_project), dual=[ ");
            for(i=0;i<n && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",x[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
    
    /* Compute norm of solution */
    norm = LPnorm(x,n,p);

    #ifdef DEBUG
        fprintf(DEBUG_FILE,"(LPp_project) p=%lf norm=%lf\n",p,norm);
    #endif
    
    /* Free and return */
    FREE
    return 1;
    
    #undef FREE
    #undef CANCEL
}

/**
    Solves the linear problem with Lp ball constraint

        argmin_s  s · z  s.t.  ||s||_p <= lambda

    This is done in O(n) time.
    
    @argument z linear weights of the problem.
    @argument n length of vector z.
    @argument p norm of the Lp ball constraint.
    @argument lambda radius of the Lp ball constraint.
    @argument s pointer to array in which to store the solution.
*/
void solveLinearLP(double *z, int n, double p, double lambda, double *s) {
    int i;
    
    // Special case where the norm p is very large
    if ( p >= LPPROJ_PLARGE ) {
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(solveLinearLP) Using large p approximation for p=%lg\n",p);
        #endif 
        // The solution approximately lies at the opposite corner of the linf ball
        for ( i = 0 ; i < n ; i++ )
            s[i] = - lambda * sign(z[i]);
        return;
    // Special case where the norm p is very close to 1
    } else if ( p <= LPPROJ_PSMALL ) {
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(solveLinearLP) Using small p approximation for p=%lg\n",p);
        #endif 
        // The solution approximately lies at the opposite corner of the l1 ball
        double largest = 0, val;
        int ilargest;
        
        // Find largest entry of z, in absolute value
        for ( i = 0 ; i < n ; i++ ) {
            val = fabs(z[i]);
            if ( val > largest) {
                largest = val;
                ilargest = i;
            }
        }
        
        // Initialize solution to all zeros
        memset(s , 0 , sizeof(double)*n);
        // Set largest entry to corner value
        s[ilargest] = - lambda * sign(z[ilargest]);
        return;
    }
    
    // General case
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"(solveLinearLP) Solving general case for p=%lg\n",p);
    #endif 

    // Compute dual norm q
    double q = 1./(1. - 1./p);
    double qMinus = q - 1;
    
    // Compute value of the dual norm for z
    double zNorm = LPnorm( z , n , q );
    
    // Get direction of solution (z is normalized to prevent numerical errors for large values of q)
    for ( i = 0 ; i < n ; i++ )
        s[i] = - sign(z[i]) * pow( fabs ( z[i] / zNorm ) , qMinus);
        
    // Compute primal norm of resultant vector
    double sNorm = LPnorm( s , n , p );
        
    // Scale back to Lp ball
    for ( i = 0 ; i < n ; i++ )
        s[i] = s[i] / sNorm * lambda;
}

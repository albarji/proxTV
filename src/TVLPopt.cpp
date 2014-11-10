/**
    Optimizers for problems dealing with general TV-Lp norm regularization.
    
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

/*  GP_TVp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_p .
        
    To do so a Gradient Projection method is used to solve its dual problem, employing
    an Lp proximity solver based on Projected Newton.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the TV norm.
*/
int GP_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws){
    double *w=NULL,*aux=NULL,*aux2=NULL,*g;
    double q,tmp,stop,dual,bestdual,lambdaMax,lambdaIni,lambdaCurrent,mu,musqrt,beta;
    int iter,stuck,nn,i,lambdaStep;
    Workspace *wsinner=NULL;
    lapack_int one=1,rc,nnp;
    
    /* Problem constants */
    #define L 4 // Lipschitz constant
    #define Linv 0.25 // Inverse of Lipschitz constant
    #define Lsqrt 2 // Squared root of Lipschitz constant
    
    #define FREE \
        if(!ws){ \
            if(w) free(w); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
        } \
        freeWorkspace(wsinner);
        
    #define CANCEL(txt,info) \
        printf("GP_TVp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Gradient to dual gap, considering also the special cases p~=1 and p~=inf */
    #define GRAD2GAP(w,g,gap,lambda,p,n,i,tmp) \
        gap = lambda * LPnorm(g, n, p); \
        for ( i = 0 ; i < n ; i++ ) \
            gap += w[i] * g[i]; \
        gap = fabs(gap);
    
    nn = n-1;
    /* Alloc memory if no workspace available */
    if(!ws){
        w = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
        aux2 = (double*)malloc(sizeof(double)*nn);
    }
    /* If a workspace is available, assign pointers */
    else{
        resetWorkspace(ws);
        w = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
    }
    if(!w || !aux || !aux2)
        {CANCEL("out of memory",info)}
    /* Gradient will be stored in an auxiliar vector */
    g = aux2;
    
    /* Alloc an additional workspace for the inner Lp ball projection solver */
    wsinner = newWorkspace(n);
    if(!wsinner)
        {CANCEL("out of memory",info)}
    
    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
        
    #ifdef TIMING
        clock_t tIni = clock();
        DUAL2PRIMAL(w,x,i)
        PRIMAL2GRAD(x,g,i)
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        fprintf(DEBUG_FILE,"%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
    #endif
    
    /* Compute dual norm */
    q = 1./(1. - 1./p);
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    nnp=nn;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = LPnorm(w, nn, q);

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax);
    #endif
    /* If it is, solution is obtained in closed form */
    if(lambda >= lambdaMax){
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* Initialize lambda step */
    if(LAMBDA_STEPS_TVLP > 1) lambdaCurrent = lambdaIni = lambda/LAMBDA_REDUCTION_TVLP;
    else lambdaCurrent = lambdaIni = lambda;
    lambdaStep = 0;
    
    /* Project point to feasible region (in case it was not feasible) */
    resetWorkspace(wsinner);
    if(!LPp_project(w,lambdaCurrent,aux,info,nn,q,wsinner))
        {CANCEL("error when invoking Lp ball projection subroutine",info)}
    for(i=0;i<nn;i++) w[i] = aux[i];
    
    /* Compute mu-convexity */
    mu = 2. - 2. * cos ( M_PI / (nn+1));
    musqrt = sqrt(mu);
    /* Compute beta */
    beta = (Lsqrt - musqrt) / (Lsqrt + musqrt);
    
    /* Start Gradient Projections iterations */
    stop = DBL_MAX; bestdual = DBL_MAX; stuck = 0;
    /* Stop when one these conditions are met:
        - Maximum number of iterations reached.
        - Dual gap small.
        - Dual gap not improving for a number of iterations. */
    for(iter=1 ; iter < MAX_ITERS_TVLP && stop > STOP_TVLP && stuck < MAX_NOIMP_TVLP ; iter++){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) Iter %d, w=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
    
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) Iter %d, g=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",g[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Lipschitz step */
        for(i=0;i<nn;i++)
            aux[i] = w[i] - Linv * g[i];
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) Iter %d, xs=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
            
        /* Projection onto lp-ball */
        resetWorkspace(wsinner);
        if(!PN_LPp(aux,lambdaCurrent,w,info,nn,p,wsinner,0,OBJGAP_LPPROX_TVLP))
            {CANCEL("error when invoking Lp ball projection subroutine",info)}

        for(i=0;i<nn;i++) w[i] = aux[i] - w[i];
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) Iter %d, xsproy=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)    
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp)   
        
        /* Check improvement in dual */
        DUALVAL(w,y,dual,i)
        if(dual < bestdual){
            bestdual = dual;
            stuck = 0;
        }
        else
            stuck++;
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) iter %d lamStep=%d lam=%lf gap=%lf, dual=%lf, bestdual=%lf\n",iter,lambdaStep,lambdaCurrent,stop,dual,bestdual); fflush(DEBUG_FILE);
        #endif
        
        /* Check dual gap */
        while(stop < STOP_TVLP){
            /* If met, move to following lambda step, if any */
            lambdaStep++;
            if(lambdaStep < LAMBDA_STEPS_TVLP){
                lambdaCurrent = pow(10, log10(lambdaIni) + lambdaStep * log10(LAMBDA_REDUCTION_TVLP) / ((double)(LAMBDA_STEPS_TVLP-1)) );
                GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp) 
                stuck = 0; bestdual = DBL_MAX;
            }
            else break;
        }
        
        #ifdef TIMING
            printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
        #endif
    }
    
    /* Termination check */
    if(iter >= MAX_ITERS_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_TVLP);
        #endif
        if(info)
            info[INFO_RC] = RC_ITERS;
    }
    else if(stuck >= MAX_NOIMP_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GP_TVp) WARNING: search stuck, improvement is not possible.\n");
        #endif
        if(info)
            info[INFO_RC] = RC_STUCK;
    }
    else if(info) info[INFO_RC] = RC_OK;

    /* Store statistics */
    if(info){
        info[INFO_ITERS] = iter-1;
        info[INFO_GAP] = fabs(stop);
    }
    
    /* Free resources and return */
    FREE
    return 1;   
    
    #undef L
    #undef Linv
    #undef Lsqrt
    #undef FREE
    #undef CANCEL
    #undef GRAD2GAP
}

/*  OGP_TVp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_p .
        
    To do so an Optimal Gradient Projection method is used to solve its dual problem, employing
    an Lp proximity solver based on Projected Newton.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the TV norm.
*/
int OGP_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws){
    double *w=NULL,*aux=NULL,*aux2=NULL,*g;
    double q,tmp,stop,dual,bestdual,lambdaMax,lambdaIni,lambdaCurrent,mu,musqrt,beta;
    int iter,stuck,nn,i,lambdaStep;
    Workspace *wsinner=NULL;
    lapack_int one=1,rc,nnp;
    
    /* Problem constants */
    #define L 4 // Lipschitz constant
    #define Linv 0.25 // Inverse of Lipschitz constant
    #define Lsqrt 2 // Squared root of Lipschitz constant
    
    #define FREE \
        if(!ws){ \
            if(w) free(w); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
        } \
        freeWorkspace(wsinner);
        
    #define CANCEL(txt,info) \
        printf("OGP_TVp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Gradient to dual gap, considering also the special cases p~=1 and p~=inf */
    #define GRAD2GAP(w,g,gap,lambda,p,n,i,tmp) \
        gap = tmp = 0; \
        if(p == mxGetInf() || p >= LPPROJ_PLARGE){ \
            for(i=0;i<n;i++){ \
                if(fabs(g[i]) > tmp) tmp = fabs(g[i]); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * tmp; \
        } \
        else if(p <= LPPROJ_PSMALL){ \
            for(i=0;i<n;i++){ \
                tmp += fabs(g[i]); \
                gap += w[i] * g[i]; \
            }\
            gap += lambda * tmp; \
        } \
        else{ \
            for(i=0;i<n;i++){ \
                tmp += pow(fabs(g[i]),p); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * pow(tmp,1.0/p); \
        } \
        gap = fabs(gap);
    
    nn = n-1;
    /* Alloc memory if no workspace available */
    if(!ws){
        w = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
        aux2 = (double*)malloc(sizeof(double)*nn);
    }
    /* If a workspace is available, assign pointers */
    else{
        resetWorkspace(ws);
        w = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
    }
    if(!w || !aux || !aux2)
        {CANCEL("out of memory",info)}
    /* Gradient will be stored in an auxiliar vector */
    g = aux2;
    
    /* Alloc an additional workspace for the inner Lp ball projection solver */
    wsinner = newWorkspace(n);
    if(!wsinner)
        {CANCEL("out of memory",info)}
    
    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
        
    #ifdef TIMING
        clock_t tIni = clock();
        DUAL2PRIMAL(w,x,i)
        PRIMAL2GRAD(x,g,i)
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
    #endif
    
    /* Compute dual norm */
    q = 1./(1. - 1./p);
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    nnp=nn;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = LPnorm(w, nn, q);

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax); fflush(DEBUG_FILE);
    #endif
    /* If it is, solution is obtained in closed form */
    if(lambda >= lambdaMax){
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* Initialize lambda step */
    if(LAMBDA_STEPS_TVLP > 1) lambdaCurrent = lambdaIni = lambda/LAMBDA_REDUCTION_TVLP;
    else lambdaCurrent = lambdaIni = lambda;
    lambdaStep = 0;
    
    /* Project point to feasible region (in case it was not feasible) */
    resetWorkspace(wsinner);
    if(!LPp_project(w,lambdaCurrent,aux,info,nn,q,wsinner))
        {CANCEL("error when invoking Lp ball projection subroutine",info)}
    for(i=0;i<nn;i++) w[i] = aux[i];
    
    /* Alloc auxiliary memory */
    double *z = (double*)malloc(sizeof(double)*nn);
    for (i=0;i<nn;i++) z[i] = w[i];
    
    /* Compute mu-convexity */
    mu = 2. - 2. * cos ( M_PI / (nn+1));
    musqrt = sqrt(mu);
    /* Compute beta */
    beta = (Lsqrt - musqrt) / (Lsqrt + musqrt);
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"mu = %lg, beta = %lg\n",mu, beta); fflush(DEBUG_FILE);
    #endif
    
    /* Start Gradient Projections iterations */
    stop = DBL_MAX; bestdual = DBL_MAX; stuck = 0;
    /* Stop when one these conditions are met:
        - Maximum number of iterations reached.
        - Dual gap small.
        - Dual gap not improving for a number of iterations. */
    for(iter=1 ; iter < MAX_ITERS_TVLP && stop > STOP_TVLP && stuck < MAX_NOIMP_TVLP ; iter++){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, w=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
    
        /* Compute the primal solution */
        DUAL2PRIMAL(z,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, g=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",g[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
        
        /* Lipschitz step */
        for(i=0;i<nn;i++)
            aux[i] = z[i] - Linv * g[i];
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, xs=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
            
        /* Projection onto lp-ball */
        resetWorkspace(wsinner);
        if(!PN_LPp(aux,lambdaCurrent,z,info,nn,p,wsinner,0,OBJGAP_LPPROX_TVLP))
            {CANCEL("error when invoking Lp ball projection subroutine",info)}
        for(i=0;i<nn;i++) z[i] = aux[i] - z[i];
            
        /* mu-convexity step */
        for(i=0;i<nn;i++) {
            // Take projected value
            tmp = z[i];
            // Compute new z value (mu-convexity step)
            z[i] = tmp + beta * (tmp - w[i]);
            // Store projected value for next update
            w[i] = tmp;
        }
        
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)    
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp)   
        
        /* Check improvement in dual */
        DUALVAL(w,y,dual,i)
        if(dual < bestdual){
            bestdual = dual;
            stuck = 0;
        }
        else
            stuck++;
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(OGP_TVp) iter=%d lamStep=%d lam=%lf gap=%lf dual=%lf bestdual=%lf\n",iter,lambdaStep,lambdaCurrent,stop,dual,bestdual); fflush(DEBUG_FILE);
        #endif
        
        /* Check dual gap */
        while(stop < STOP_TVLP){
            /* If met, move to following lambda step, if any */
            lambdaStep++;
            if(lambdaStep < LAMBDA_STEPS_TVLP){
                lambdaCurrent = pow(10, log10(lambdaIni) + lambdaStep * log10(LAMBDA_REDUCTION_TVLP) / ((double)(LAMBDA_STEPS_TVLP-1)) );
                GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp) 
                stuck = 0; bestdual = DBL_MAX;
            }
            else break;
        }
        
        #ifdef TIMING
            printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
        #endif
    }
    
    /* Termination check */
    if(iter >= MAX_ITERS_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(OGP_TVp) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_TVLP); fflush(DEBUG_FILE);
        #endif
        if(info)
            info[INFO_RC] = RC_ITERS;
    }
    else if(stuck >= MAX_NOIMP_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(OGP_TVp) WARNING: search stuck, improvement is not possible.\n"); fflush(DEBUG_FILE);
        #endif
        if(info)
            info[INFO_RC] = RC_STUCK;
    }
    else if(info) info[INFO_RC] = RC_OK;

    /* Store statistics */
    if(info){
        info[INFO_ITERS] = iter-1;
        info[INFO_GAP] = fabs(stop);
    }
    
    /* Free resources and return */
    FREE
    free(z); // Free auxiliary memory
    return 1;   
    
    #undef L
    #undef Linv
    #undef Lsqrt
    #undef FREE
    #undef CANCEL
    #undef GRAD2GAP
}

/*  FISTA_TVp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_p .
        
    To do so a FISTA method is used to solve its dual problem, employing
    an Lp proximity solver based on Projected Newton.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the TV norm.
*/
int FISTA_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws){
    double *w=NULL,*aux=NULL,*aux2=NULL,*g;
    double q,tmp,stop,dual,bestdual,lambdaMax,lambdaIni,lambdaCurrent,mu,musqrt,beta;
    int iter,stuck,nn,i,lambdaStep;
    Workspace *wsinner=NULL;
    lapack_int one=1,rc,nnp;
    
    /* Problem constants */
    #define L 4 // Lipschitz constant
    #define Linv 0.25 // Inverse of Lipschitz constant
    #define Lsqrt 2 // Squared root of Lipschitz constant
    
    #define FREE \
        if(!ws){ \
            if(w) free(w); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
        } \
        freeWorkspace(wsinner);
        
    #define CANCEL(txt,info) \
        printf("FISTA_TVp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Gradient to dual gap, considering also the special cases p~=1 and p~=inf */
    #define GRAD2GAP(w,g,gap,lambda,p,n,i,tmp) \
        gap = tmp = 0; \
        if(p == mxGetInf() || p >= LPPROJ_PLARGE){ \
            for(i=0;i<n;i++){ \
                if(fabs(g[i]) > tmp) tmp = fabs(g[i]); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * tmp; \
        } \
        else if(p <= LPPROJ_PSMALL){ \
            for(i=0;i<n;i++){ \
                tmp += fabs(g[i]); \
                gap += w[i] * g[i]; \
            }\
            gap += lambda * tmp; \
        } \
        else{ \
            for(i=0;i<n;i++){ \
                tmp += pow(fabs(g[i]),p); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * pow(tmp,1.0/p); \
        } \
        gap = fabs(gap);
    
    nn = n-1;
    /* Alloc memory if no workspace available */
    if(!ws){
        w = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
        aux2 = (double*)malloc(sizeof(double)*nn);
    }
    /* If a workspace is available, assign pointers */
    else{
        resetWorkspace(ws);
        w = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
    }
    if(!w || !aux || !aux2)
        {CANCEL("out of memory",info)}
    /* Gradient will be stored in an auxiliar vector */
    g = aux2;
    
    /* Alloc an additional workspace for the inner Lp ball projection solver */
    wsinner = newWorkspace(n);
    if(!wsinner)
        {CANCEL("out of memory",info)}
    
    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
        
    #ifdef TIMING
        clock_t tIni = clock();
        DUAL2PRIMAL(w,x,i)
        PRIMAL2GRAD(x,g,i)
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
    #endif
    
    /* Compute dual norm */
    q = 1./(1. - 1./p);
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    nnp=nn;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = LPnorm(w, nn, q);

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax); fflush(DEBUG_FILE);
    #endif
    /* If it is, solution is obtained in closed form */
    if(lambda >= lambdaMax){
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* Initialize lambda step */
    if(LAMBDA_STEPS_TVLP > 1) lambdaCurrent = lambdaIni = lambda/LAMBDA_REDUCTION_TVLP;
    else lambdaCurrent = lambdaIni = lambda;
    lambdaStep = 0;
    
    /* Project point to feasible region (in case it was not feasible) */
    resetWorkspace(wsinner);
    if(!LPp_project(w,lambdaCurrent,aux,info,nn,q,wsinner))
        {CANCEL("error when invoking Lp ball projection subroutine",info)}
    for(i=0;i<nn;i++) w[i] = aux[i];
    
    /* Alloc auxiliary memory */
    double *z = (double*)malloc(sizeof(double)*nn);
    double *wpre = (double*)malloc(sizeof(double)*nn);
    for (i=0;i<nn;i++) z[i] = w[i];

    double fistaStep = 1, fistaStepPrev;
    
    /* Start Gradient Projections iterations */
    stop = DBL_MAX; bestdual = DBL_MAX; stuck = 0;
    /* Stop when one these conditions are met:
        - Maximum number of iterations reached.
        - Dual gap small.
        - Dual gap not improving for a number of iterations. */
    for(iter=1 ; iter < MAX_ITERS_TVLP && stop > STOP_TVLP && stuck < MAX_NOIMP_TVLP ; iter++){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, w=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
    
        /* Compute the primal solution */
        DUAL2PRIMAL(z,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, g=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",g[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
        
        /* Lipschitz step */
        for(i=0;i<nn;i++)
            aux[i] = z[i] - Linv * g[i];
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"Iter %d, xs=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
            fprintf(DEBUG_FILE,"]\n"); fflush(DEBUG_FILE);
        #endif
        
        /* Store current point */
        for(i=0;i<nn;i++) {
            wpre[i] = w[i];
        }
        
        /* Projection onto lq-ball */
        resetWorkspace(wsinner);
        if(!PN_LPp(aux,lambdaCurrent,w,info,nn,p,wsinner,0,OBJGAP_LPPROX_TVLP))
            {CANCEL("error when invoking Lp ball projection subroutine",info)}
        for(i=0;i<nn;i++) w[i] = aux[i] - w[i];
            
        /* FISTA readjustment step */
        fistaStepPrev = fistaStep;
        fistaStep = (1. + sqrt(1.+4*fistaStep*fistaStep)) / 2.;
        beta = (fistaStepPrev-1.) / fistaStep;
        for(i=0;i<nn;i++) {
            // Compute new z value (FISTA step)
            z[i] = w[i] + beta * (w[i] - wpre[i]);
        }
        
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)    
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp)   
        
        /* Check improvement in dual */
        DUALVAL(w,y,dual,i)
        if(dual < bestdual){
            bestdual = dual;
            stuck = 0;
        }
        else
            stuck++;
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FISTA_TVp) iter=%d lamStep=%d lam=%lf gap=%lf dual=%lf bestdual=%lf\n",iter,lambdaStep,lambdaCurrent,stop,dual,bestdual); fflush(DEBUG_FILE);
        #endif
        
        /* Check dual gap */
        while(stop < STOP_TVLP){
            /* If met, move to following lambda step, if any */
            lambdaStep++;
            if(lambdaStep < LAMBDA_STEPS_TVLP){
                lambdaCurrent = pow(10, log10(lambdaIni) + lambdaStep * log10(LAMBDA_REDUCTION_TVLP) / ((double)(LAMBDA_STEPS_TVLP-1)) );
                GRAD2GAP(w,g,stop,lambdaCurrent,p,n,i,tmp) 
                stuck = 0; bestdual = DBL_MAX;
            }
            else break;
        }
        
        #ifdef TIMING
            printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
        #endif
    }
    
    /* Termination check */
    if(iter >= MAX_ITERS_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FISTA_TVp) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_TVLP); fflush(DEBUG_FILE);
        #endif
        if(info)
            info[INFO_RC] = RC_ITERS;
    }
    else if(stuck >= MAX_NOIMP_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FISTA_TVp) WARNING: search stuck, improvement is not possible.\n"); fflush(DEBUG_FILE);
        #endif
        if(info)
            info[INFO_RC] = RC_STUCK;
    }
    else if(info) info[INFO_RC] = RC_OK;

    /* Store statistics */
    if(info){
        info[INFO_ITERS] = iter-1;
        info[INFO_GAP] = fabs(stop);
    }
    
    /* Free resources and return */
    FREE
    free(z); free(wpre); // Free auxiliary memory
    return 1;   
    
    #undef L
    #undef Linv
    #undef Lsqrt
    #undef FREE
    #undef CANCEL
    #undef GRAD2GAP
}

/*  FW_TVp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_p .
        
    To do so a Franke-Wolfe conditional gradient is used to solve its dual problem
    
        min_u   1/2 ||D^T u||^2_2 - u^T D y    s.t. ||u||_q <= lambda
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the TV norm.
*/
int FW_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws){
    int nn, iter, i;
    double *w=NULL, *aux=NULL, *aux2=NULL, *g=NULL;
    double stop, q, lambdaMax, gap, step, den, tmp, fval, fvalPrev;
    lapack_int one=1,rc,nnp;
    
    /* Gradient to dual gap, considering also the special cases p~=1 and p~=inf */
    #define GRAD2GAP(w,g,gap,lambda,p,n,i,tmp) \
        gap = tmp = 0; \
        if(p == mxGetInf() || p >= LPPROJ_PLARGE){ \
            for(i=0;i<n;i++){ \
                if(fabs(g[i]) > tmp) tmp = fabs(g[i]); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * tmp; \
        } \
        else if(p <= LPPROJ_PSMALL){ \
            for(i=0;i<n;i++){ \
                tmp += fabs(g[i]); \
                gap += w[i] * g[i]; \
            }\
            gap += lambda * tmp; \
        } \
        else{ \
            for(i=0;i<n;i++){ \
                tmp += pow(fabs(g[i]),p); \
                gap += w[i] * g[i]; \
            } \
            gap += lambda * pow(tmp,1.0/p); \
        } \
        gap = fabs(gap);
        
    #define FREE \
        if (!ws) { \
            if(w) free(w); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
            if(g) free(g); \
        }
        
    #define CANCEL(txt,info) \
        printf("FW_TVp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;

    /* Prepare memory/workspace */
    nn = n-1;
    if ( ws ) {
        w = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
        g = getDoubleWorkspace(ws);
    } else {
        w = (double*) malloc( sizeof(double) * nn );
        aux = (double*) malloc( sizeof(double) * nn );
        aux2 = (double*) malloc( sizeof(double) * nn );
        g = (double*) malloc( sizeof(double) * nn );
    }
    if (!w || !aux || !aux2 || !g)
        {CANCEL("Error retrieving memory from workspace",info)}

    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
    
    /* Compute dual norm */
    q = 1./(1. - 1./p);
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    nnp=nn;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = LPnorm(w, nn, q);

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax);
    #endif
    /* If it is, solution is obtained in closed form */
    if(lambda >= lambdaMax){
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambda,p,nn,i,tmp)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* Else, use radially projected point as starting guess */
    for(i=0;i<nn;i++)
        w[i] = w[i] / lambdaMax * lambda;
        
    /* Initial dual value */
    DUALVAL(w,y,fval,i);
    fvalPrev = DBL_MAX;

    /* Franke-Wolfe loop */
    gap = step = DBL_MAX;
    iter = 0;
    while(gap > STOP_TVLP && (fvalPrev - fval) > MIN_IMP_TVLP && iter < MAX_ITERS_TVLPFW) {
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, w=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
    
        /* Get gradient */
        DUAL2PRIMAL(w,x,i)
        PRIMAL2GRAD(x,g,i)
        
        /* Store current primal value */
        fvalPrev = fval;
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, g=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",g[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
    
        /* Solve linear surrogate problem: argmin_s  s · gradient(w)  s.t.  ||s||_q <= lambda */
        solveLinearLP(g,nn,q,lambda,aux);
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, s=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Get displacement vector */
        for ( i = 0 ; i < nn ; i++ )
            aux2[i] = aux[i] - w[i];
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, disp=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux2[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Compute surrogate dual gap */
        gap = 0;
        for ( i = 0 ; i < nn ; i++ )
            gap -= aux2[i] * g[i];
        
        /* Compute unbounded minimizing stepsize: gap / ||D' · (w-s)|| */
        /* (D'd)[0] = -d[0], (D'd)[i] = sum(d[i-1] - d[i]), (D'd)[nn] = d[nn] */
        den = aux2[0] * aux2[0];
        for ( i = 0 ; i < nn-1 ; i++ ) {
            tmp = aux2[i] - aux2[i+1];
            den += tmp * tmp;
        }
        den += aux2[nn] * aux2[nn];
        step = gap / den;
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, step=%lf\n",iter,step);
        #endif
        
        /* Clip down to interval [0,1] */
        step = min(max(0.,step),1.);
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) Iter %d, stepClip=%lf\n",iter,step);
        #endif
        
        /* Peform a step along the surrogate solution */
        for ( i = 0 ; i < nn ; i++ )
            w[i] = (1.-step) * w[i] + step * aux[i];
            
        /* Compute new dual value */
        DUALVAL(w,y,fval,i);
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(FW_TVp) surrgap=%lf, dual=%lf\n",gap,fval); fflush(stdout);
        #endif
            
        /* If surrogate gap is met, check whether the same applies for the real gap */
        if ( gap <= STOP_TVLP ) {
            /* Compute the primal solution */
            DUAL2PRIMAL(w,x,i)
            /* Gradient evaluation */
            PRIMAL2GRAD(x,g,i)
            /* Compute real dual gap */
            GRAD2GAP(w,g,gap,lambda,p,nn,i,tmp)
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(FW_TVp) gap=%lf\n",gap); fflush(stdout);
            #endif
        }
        
        iter++;
    }
        
    /* Return info if requested */
    if (info) {
        info[INFO_RC] = RC_OK;
        info[INFO_ITERS] = iter;
        info[INFO_GAP] = gap;
    }
    
    #undef FREE
    #undef CANCEL
    #undef GRAD2GAP
}

/*  GPFW_TVp

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_p .
        
    To do so a hybrid Gradient Projection + Frank-Wolfe method is used to solve its dual problem, employing
    an Lp proximity solver based on Projected Newton. Iterations of both methods are interleaved for
    further efficiency.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - p: degree of the TV norm.
*/
int GPFW_TVp(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws) {
    double *w=NULL,*aux=NULL,*aux2=NULL,*g=NULL;
    double q,tmp,stop,dual,bestdual,lambdaMax,num,den,step;
    int iter,stuck,nn,i,lambdaStep,cycle;
    Workspace *wsinner=NULL;
    lapack_int one=1,rc,nnp;
    
    /* Problem constants */
    #define Linv 0.25 // Inverse of Lipschitz constant
    
    #define FREE \
        if(!ws){ \
            if(w) free(w); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
            if(g) free(g); \
        } \
        freeWorkspace(wsinner);
        
    #define CANCEL(txt,info) \
        printf("GPFW_TVp: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
        
    /* Gradient to dual gap, considering also the special cases p~=1 and p~=inf */
    #define GRAD2GAP(w,g,gap,lambda,p,n,i) \
        gap = lambda * LPnorm(g, n, p); \
        for ( i = 0 ; i < n ; i++ ) \
            gap += w[i] * g[i]; \
        gap = fabs(gap);
        
    /* If p > 10 then FW is terrible: switch to pure GP */
    if ( p > 10 )
        return GP_TVp(y, lambda, x, info, n, p, ws);
    
    nn = n-1;
    /* Alloc memory if no workspace available */
    if(!ws){
        w = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
        aux2 = (double*)malloc(sizeof(double)*nn);
        g = (double*)malloc(sizeof(double)*nn);
    }
    /* If a workspace is available, assign pointers */
    else{
        resetWorkspace(ws);
        w = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
        g = getDoubleWorkspace(ws);
    }
    if(!w || !aux || !aux2 || !g)
        {CANCEL("out of memory",info)}
    /* Gradient will be stored in an auxiliar vector */
    
    /* Alloc an additional workspace for the inner Lp ball projection solver */
    wsinner = newWorkspace(n);
    if(!wsinner)
        {CANCEL("out of memory",info)}
    
    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
        
    #ifdef TIMING
        clock_t tIni = clock();
        DUAL2PRIMAL(w,x,i)
        PRIMAL2GRAD(x,g,i)
        GRAD2GAP(w,g,stop,lambda,p,nn,i)
        fprintf(DEBUG_FILE,"%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
    #endif
    
    /* Compute dual norm */
    q = 1./(1. - 1./p);
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    nnp=nn;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = LPnorm(w, nn, q);

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"(GPFW_TVp) lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax);
    #endif
    /* If it is, solution is obtained in closed form */
    if(lambda >= lambdaMax){
        /* Compute the primal solution */
        DUAL2PRIMAL(w,x,i)
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(w,g,stop,lambda,p,nn,i)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* Project point to feasible region (in case it was not feasible) */
    resetWorkspace(wsinner);
    if(!LPp_project(w,lambda,aux,info,nn,q,wsinner))
        {CANCEL("error when invoking Lp ball projection subroutine",info)}
    for(i=0;i<nn;i++) w[i] = aux[i];
    
    /* Compute initial primal vector */
    DUAL2PRIMAL(w,x,i)
    /* Initial gradient evaluation */
    PRIMAL2GRAD(x,g,i)
    
    /* Start Gradient Projections iterations */
    stop = DBL_MAX; bestdual = DBL_MAX; stuck = 0; cycle = 0;
    /* Stop when one these conditions are met:
        - Maximum number of iterations reached.
        - Dual gap small.
        - Dual gap not improving for a number of iterations. */
    for(iter=1 ; iter < MAX_ITERS_TVLPGPFW && stop > STOP_TVLP && stuck < MAX_NOIMP_TVLP_GPFW ; iter++, cycle++){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, w=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",w[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, g=[ ",iter);
            for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",g[i]);
            fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /*** Gradient Projection cycle ***/
        if ( cycle < 0 || !(cycle % FW_CYCLES_TVLP) ) {
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) GP iteration\n");
            #endif
                
            /* Lipschitz step */
            for(i=0;i<nn;i++)
                aux2[i] = w[i] - Linv * g[i];
                
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, xs=[ ",iter);
                for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux2[i]);
                fprintf(DEBUG_FILE,"]\n");
            #endif
                
            /* Projection onto lp-ball */
            resetWorkspace(wsinner);
            if(!PN_LPp(aux2,lambda,aux,info,nn,p,wsinner,0,OBJGAP_LPPROX_TVLP))
                {CANCEL("error when invoking Lp ball projection subroutine",info)}
            for(i=0;i<nn;i++) aux[i] = aux2[i] - aux[i];
                
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, xsproy=[ ",iter);
                for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
                fprintf(DEBUG_FILE,"]\n");
            #endif
            
            /* Find optimal stepsize between original w and Lipschitz+Projected point */
            
            /* Get displacement vector */
            for ( i = 0 ; i < nn ; i++ )
                aux2[i] = aux[i] - w[i];
            
            /* Compute stepsize numerator */
            num = 0;
            for ( i = 0 ; i < nn ; i++ )
                num -= aux2[i] * g[i];
            
            /* Compute unbounded minimizing stepsize: num / ||D' · d|| */
            /* (D'd)[0] = -d[0], (D'd)[i] = sum(d[i-1] - d[i]), (D'd)[nn] = d[nn] */
            den = aux2[0] * aux2[0];
            for ( i = 0 ; i < nn-1 ; i++ ) {
                tmp = aux2[i] - aux2[i+1];
                den += tmp * tmp;
            }
            den += aux2[nn] * aux2[nn];
            step = num / den;
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) step=%lf\n",step); fflush(stdout);
            #endif
            
            /* If we get negative stepsize, don't clip it. We will enworse the obj. function in this step,
               but this should help moving out of numericaly unstable regions */
            if ( step < 0 || cycle < 0 ) {
                step = 1;
                cycle = -MAX_ITERS_TVLPFW;
            } else            
                /* Clip down to interval [0,1] */
            step = min(max(0.,step),1.);
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) stepClip=%lf\n",step); fflush(stdout);
            #endif
            
            /* Peform a step along the surrogate solution */
            if ( step == 1 )
                memcpy(w, aux, sizeof(double)*nn);
            else
                for ( i = 0 ; i < nn ; i++ )
                    w[i] = (1.-step) * w[i] + step * aux[i];

            /* Compute the primal solution */
            DUAL2PRIMAL(w,x,i)
            /* Gradient evaluation */
            PRIMAL2GRAD(x,g,i)
            /* Compute real dual gap */
            GRAD2GAP(w,g,stop,lambda,p,nn,i)
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) gap=%lf\n",stop); fflush(stdout);
            #endif
        
        /*** Frank-Wolfe step ***/
        } else {
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) FW iteration\n");
            #endif
    
            /* Solve linear surrogate problem: argmin_s  s · gradient(w)  s.t.  ||s||_q <= lambda */
            solveLinearLP(g,nn,q,lambda,aux);
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, s=[ ",iter);
                for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux[i]);
                fprintf(DEBUG_FILE,"]\n");
            #endif
            
            /* Get displacement vector */
            for ( i = 0 ; i < nn ; i++ )
                aux2[i] = aux[i] - w[i];
                
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) Iter %d, disp=[ ",iter);
                for(i=0;i<nn && i<DEBUG_N;i++) fprintf(DEBUG_FILE,"%g ",aux2[i]);
                fprintf(DEBUG_FILE,"]\n");
            #endif
            
            /* Compute surrogate dual gap */
            stop = 0;
            for ( i = 0 ; i < nn ; i++ )
                stop -= aux2[i] * g[i];
            
            /* Compute unbounded minimizing stepsize: gap / ||D' · (w-s)|| */
            /* (D'd)[0] = -d[0], (D'd)[i] = sum(d[i-1] - d[i]), (D'd)[nn] = d[nn] */
            den = aux2[0] * aux2[0];
            for ( i = 0 ; i < nn-1 ; i++ ) {
                tmp = aux2[i] - aux2[i+1];
                den += tmp * tmp;
            }
            den += aux2[nn] * aux2[nn];
            step = stop / den;
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) step=%lf\n",step); fflush(stdout);
            #endif
            
            /* Clip down to interval [0,1] */
            step = min(max(0.,step),1.);
            
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"(GPFW_TVp) stepClip=%lf\n",step); fflush(stdout);
            #endif
            
            /* Peform a step along the surrogate solution */
            for ( i = 0 ; i < nn ; i++ )
                w[i] = (1.-step) * w[i] + step * aux[i];
            
            /* Compute the primal solution */
            DUAL2PRIMAL(w,x,i)
            /* Gradient evaluation */
            PRIMAL2GRAD(x,g,i)
                
            /* If surrogate gap is met, check whether the same applies for the real gap */
            if ( stop <= STOP_TVLP ) {
                /* Compute real dual gap */
                GRAD2GAP(w,g,stop,lambda,p,nn,i)
            } 
        }
        
        /* Check improvement in dual */
        DUALVAL(w,y,dual,i)
        /* Check significant improvement in dual */
        if(dual < bestdual - MIN_IMP_TVLP){
            stuck = 0;
            /* If FW cycles have been banned, restore them again */
            if ( cycle < 0 ) cycle = 0;
        }
        else
            stuck++;

        if ( dual < bestdual ) 
            bestdual = dual;
            
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GPFW_TVp) iter %d gap=%lf, dual=%lf, bestdual=%lf\n",iter,stop,dual,bestdual); fflush(DEBUG_FILE); fflush(stdout);
        #endif
        
        #ifdef TIMING
            printf("%lf\t%lf\n",((double)(clock()-tIni))/((double)CLOCKS_PER_SEC),stop);
        #endif
    }
    
    /* Termination check */
    if(iter >= MAX_ITERS_TVLP){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GPFW_TVp) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_TVLP); fflush(stdout);
        #endif
        if(info)
            info[INFO_RC] = RC_ITERS;
    }
    else if(stuck >= MAX_NOIMP_TVLP_GPFW){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(GPFW_TVp) WARNING: search stuck, improvement is not possible.\n"); fflush(stdout);
        #endif
        if(info)
            info[INFO_RC] = RC_STUCK;
    }
    else if(info) info[INFO_RC] = RC_OK;

    /* Store statistics */
    if(info){
        info[INFO_ITERS] = iter-1;
        info[INFO_GAP] = fabs(stop);
    }
    
    /* Free resources and return */
    FREE
    return 1;   
    
    #undef L
    #undef FREE
    #undef CANCEL
    #undef GRAD2GAP
}


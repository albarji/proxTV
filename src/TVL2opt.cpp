/**
    Optimizers for problems dealing with TV-L2 norm regularization.
    
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

/*  more_TV2

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_2 .
        
    To do so a More-Sorensen algorithm is used to solve its dual problem.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
*/
int more_TV2(double *y,double lambda,double *x,double *info,int n){
    int nn=n-1,i;
    double stop,tmp,tmp2,lam,pNorm,qNorm,pNormSq,dist;
    double *Dy,*alpha,*beta,*minus,*p,*aux;
    lapack_int one=1,rc,nnp=nn;
    
    /* Macros */
    
    // Solves Rx = y for lower bidiagonal R given by diag. alpha and subdiag. beta using forward substitution
    // Returns the solution overwriting y
    #define FW_SUBS(alpha,beta,y,n,i) \
        y[0] /= alpha[0]; \
        for(i=1;i<n;i++) \
            y[i] = (y[i] - beta[i-1] * y[i-1]) / alpha[i];
            
    #define GRAD2GAP(w,g,gap,lambda,n,i,tmp) \
        gap = tmp = 0; \
        for(i=0;i<n;i++){ \
            tmp += g[i]*g[i]; \
            gap += w[i] * g[i]; \
        } \
        gap += lambda * sqrt(tmp); \
        gap = fabs(gap);
    
    #define FREE \
        if(Dy) free(Dy); \
        if(minus) free(minus); \
        if(alpha) free(alpha); \
        if(beta) free(beta); \
        if(p) free(p); \
        if(aux) free(aux);
        
    #define CANCEL(txt,info) \
        printf("more_TV2: %s\n",txt); \
        FREE \
        info[INFO_RC] = RC_ERROR;\
        return 0;
    
    /* Alloc memory */
    Dy = (double*)malloc(sizeof(double)*nn);
    minus = (double*)malloc(sizeof(double)*(nn-1));
    alpha = (double*)malloc(sizeof(double)*nn);
    beta = (double*)malloc(sizeof(double)*(nn-1));
    p = (double*)malloc(sizeof(double)*nn);
    aux = (double*)malloc(sizeof(double)*nn);
    if(!Dy || !minus || !alpha || !beta || !p || !aux){CANCEL("out of memory",info)}

    /* Precomputations */
    
    for(i=0;i<nn-1;i++){
        Dy[i] = -y[i] + y[i+1];
        minus[i] = -1;
    }
    Dy[nn-1] = -y[nn-1] + y[nn];
    
    /* Iterate till convergence */
    stop = DBL_MAX;
    info[INFO_ITERS] = 0;
    lam = 0;
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"--------------- Start ---------------\n",lam);
        fprintf(DEBUG_FILE,"lam=%lf\n",lam);
    #endif
    while(stop > STOP_MS && info[INFO_ITERS] < MAX_ITERS_MS){
        /* Generate tridiagonal representation of Hessian */
        tmp = 2+lam;
        for(i=0;i<nn;i++)
            alpha[i] = tmp;    
        memcpy((void*)beta,(void*)minus,sizeof(double)*(nn-1));

        /* Compute tridiagonal factorization of Hessian */
        dpttrf_(&nnp,alpha,beta,&rc);
        
        /* Obtain p by solving Cholesky system */
        memcpy((void*)aux,(void*)Dy,sizeof(double)*nn);
        dpttrs_(&nnp, &one, alpha, beta, aux, &nnp, &rc);
        memcpy((void*)p,(void*)aux,sizeof(double)*nn);
        pNorm = 0; for(i=0;i<nn;i++) pNorm += aux[i]*aux[i];
        pNormSq = sqrt(pNorm);
        
        /* Compute Cholesky matrix */
        for(i=0;i<nn-1;i++){
            alpha[i] = tmp = sqrt(alpha[i]);
            beta[i] *= tmp;
        }
        alpha[nn-1] = sqrt(alpha[nn-1]);


        /* Obtain q by solving yet another system */
        FW_SUBS(alpha,beta,aux,nn,i)
        qNorm = 0; for(i=0;i<nn;i++) qNorm += aux[i]*aux[i];

        /* Update the constraint satisfaction parameter of the algorithm */
        lam += (pNorm / qNorm) * (pNormSq - lambda) / lambda;
        /* If negative, set to zero */
        if(lam < 0) lam = 0;
        
        /* Compute distance to boundary */
        dist = pNormSq - lambda;
        /* Check if the distance criterion is met
          If we are in the lam=0 case and the p is in the interior, it is automatically met */
        if((!lam && dist <= 0) || fabs(dist) <= STOP_MSSUB){
            /* Compute dual gap */
            DUAL2PRIMAL(p,x,i)
            PRIMAL2GRAD(x,aux,i)
            GRAD2GAP(p,aux,stop,lambda,nn,i,tmp)
            stop = fabs(stop);
        }
        //else stop = dist;
        else stop = DBL_MAX;

        info[INFO_ITERS]++;
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"--------------- End of iteration %lf ---------------\n",info[INFO_ITERS]);
            fprintf(DEBUG_FILE,"p=["); for(i=0;i<nn;i++) fprintf(DEBUG_FILE,"%lf ",p[i]); fprintf(DEBUG_FILE,"]\n");
        #endif
    }
    
    info[INFO_GAP] = stop;
    
    /* Termination check */
    if(info[INFO_ITERS] >= MAX_ITERS_MS){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(more_TV2) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_PN);
        #endif
        info[INFO_RC] = RC_ITERS;
    }
    else info[INFO_RC] = RC_OK;
    
    FREE
    return 1;
    
    #undef FW_SUBS
    #undef GRAD2GAP
    #undef FREE
    #undef CANCEL
}

/*  morePG_TV2

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_2 .
        
    To do so a More-Sorensen + Projected Gradient algorithm is used to solve its dual problem.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
*/
int morePG_TV2(double *y,double lambda,double *x,double *info,int n,Workspace *ws){
    int nn=n-1,i,iters;
    double stop,tmp,tmp2,lam,pNorm,qNorm,pNormSq,dist;
    double *Dy,*alpha,*beta,*minus,*p,*aux;
    lapack_int one=1,rc,nnp=nn;
    
    /* Macros */
    
    // Solves Rx = y for lower bidiagonal R given by diag. alpha and subdiag. beta using forward substitution
    // Returns the solution overwriting y
    #define FW_SUBS(alpha,beta,y,n,i) \
        y[0] /= alpha[0]; \
        for(i=1;i<n;i++) \
            y[i] = (y[i] - beta[i-1] * y[i-1]) / alpha[i];
            
    #define GRAD2GAP(w,g,gap,lambda,n,i,tmp) \
        gap = tmp = 0; \
        for(i=0;i<n;i++){ \
            tmp += g[i]*g[i]; \
            gap += w[i] * g[i]; \
        } \
        gap += lambda * sqrt(tmp); \
        gap = fabs(gap);
        
    #define NORM(x,n,i,tmp) \
        tmp = 0; \
        for(i=0;i<n;i++) tmp += x[i]*x[i]; tmp = sqrt(tmp);
    
    #define FREE \
        if(!ws) { \
            if(Dy) free(Dy); \
            if(minus) free(minus); \
            if(alpha) free(alpha); \
            if(beta) free(beta); \
            if(p) free(p); \
            if(aux) free(aux); \
        }
        
    #define CANCEL(txt,info) \
        printf("more_TV2: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
    
    /* Alloc memory if needed */
    if(!ws){
        p = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
    }
    else{
        p = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
    }
    if(!p || !aux)
        {CANCEL("out of memory",info)}
    
    stop = DBL_MAX;
    iters = 0;
    
    /* Rule-of-thumb to check if PG might help */
    NORM(y,n,i,tmp)
    if(tmp > lambda){
        #define STEP 0.25
        #define MAX_PG 50
        
        /* Warm restart (if possible) */
        if(ws && ws->warm){
            memcpy(p,ws->warmDual,sizeof(double)*nn);
            DUAL2PRIMAL(p,x,i)
            PRIMAL2GRAD(x,aux,i)
        }
        /* Else start at 0 */
        else{
            for(i=0;i<nn;i++){
                p[i] = 0;
                aux[i] = y[i] - y[i+1];
            }
        }
        
        /* Projected Gradient iterations */
        while(stop > STOP_MS && iters < MAX_PG){
            /* Gradient step */
            for(i=0;i<nn;i++) p[i] = p[i] - STEP * aux[i];
            
            /* Projection step */
            NORM(p,n,i,tmp)
            if(tmp > lambda){
                tmp = lambda / tmp;
                for(i=0;i<nn;i++) p[i] *= tmp;
            }
            
            DUAL2PRIMAL(p,x,i)
            PRIMAL2GRAD(x,aux,i)
            GRAD2GAP(p,aux,stop,lambda,nn,i,tmp)

            iters++;
        }
        
        /* Stop if solution is good enough */
        if(stop <= STOP_MS){
            if(info){
                info[INFO_ITERS] = iters;
                info[INFO_GAP] = fabs(stop);
                info[INFO_RC] = RC_OK;
            }
            /* Store info for warm restart */
            if(ws){
                memcpy(ws->warmDual,p,sizeof(double)*nn);
                ws->warmLambda = 0;
                ws->warm = 1;
            }
            FREE
            return 1;
        }
        
        #undef STEP
        #undef MAX_PG
    }
    
    /* Alloc more memory */
    if(!ws){
        Dy = (double*)malloc(sizeof(double)*nn);
        minus = (double*)malloc(sizeof(double)*(nn-1));
        alpha = (double*)malloc(sizeof(double)*nn);
        beta = (double*)malloc(sizeof(double)*(nn-1));
    }
    else{
        Dy = getDoubleWorkspace(ws);
        minus = getDoubleWorkspace(ws);
        alpha = getDoubleWorkspace(ws);
        beta = getDoubleWorkspace(ws);
    }
    if(!Dy || !minus || !alpha || !beta || !p || !aux)
        {CANCEL("out of memory",info)}

    /* Precomputations */
    
    for(i=0;i<nn-1;i++){
        Dy[i] = -y[i] + y[i+1];
        minus[i] = -1;
    }
    Dy[nn-1] = -y[nn-1] + y[nn];
    
    /* Warm restart if possible */
    if(ws && ws->warm){
        lam = ws->warmLambda;
    }
    else lam = 0;
    
    /* Iterate till convergence */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"--------------- Start ---------------\n",lam);
        fprintf(DEBUG_FILE,"lam=%lf\n",lam);
    #endif
    while(stop > STOP_MS && iters < MAX_ITERS_MS){
        /* Generate tridiagonal representation of Hessian */
        tmp = 2+lam;
        for(i=0;i<nn;i++)
            alpha[i] = tmp;    
        memcpy((void*)beta,(void*)minus,sizeof(double)*(nn-1));

        /* Compute tridiagonal factorization of Hessian */
        dpttrf_(&nnp,alpha,beta,&rc);
        
        /* Obtain p by solving Cholesky system */
        memcpy((void*)aux,(void*)Dy,sizeof(double)*nn);
        dpttrs_(&nnp, &one, alpha, beta, aux, &nnp, &rc);
        memcpy((void*)p,(void*)aux,sizeof(double)*nn);
        pNorm = 0; for(i=0;i<nn;i++) pNorm += aux[i]*aux[i];
        pNormSq = sqrt(pNorm);
        
        /* Compute Cholesky matrix */
        for(i=0;i<nn-1;i++){
            alpha[i] = tmp = sqrt(alpha[i]);
            beta[i] *= tmp;
        }
        alpha[nn-1] = sqrt(alpha[nn-1]);


        /* Obtain q by solving yet another system */
        FW_SUBS(alpha,beta,aux,nn,i)
        qNorm = 0; for(i=0;i<nn;i++) qNorm += aux[i]*aux[i];

        /* Update the constraint satisfaction parameter of the algorithm */
        lam += (pNorm / qNorm) * (pNormSq - lambda) / lambda;
        /* If negative, set to zero */
        if(lam < 0) lam = 0;
        
        /* Compute distance to boundary */
        dist = pNormSq - lambda;
        /* Check if the distance criterion is met
          If we are in the lam=0 case and the p is in the interior, it is automatically met */
        if((!lam && dist <= 0) || fabs(dist) <= STOP_MSSUB){
            /* Compute dual gap */
            DUAL2PRIMAL(p,x,i)
            PRIMAL2GRAD(x,aux,i)
            GRAD2GAP(p,aux,stop,lambda,nn,i,tmp)
            stop = fabs(stop);
        }
        //else stop = dist;
        else stop = DBL_MAX;

        iters++;
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"stop=%lf, dist=%lf\n",stop,dist);
            fprintf(DEBUG_FILE,"--------------- End of iteration %d ---------------\n",iters);
        #endif
    }
    
    if(info){
        info[INFO_GAP] = stop;
        info[INFO_ITERS] = iters;
    }
    
    /* Termination check */
    if(iters >= MAX_ITERS_MS){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(more_TV2) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_PN);
        #endif
        if(info) info[INFO_RC] = RC_ITERS;
    }
    else if(info) info[INFO_RC] = RC_OK;
    
    /* Store info for warm restart */
    if(ws){
        memset(ws->warmDual,0,sizeof(double)*nn);
        ws->warmLambda = lam;
        ws->warm = 1;
    }
    
    FREE
    return 1;
    
    #undef FW_SUBS
    #undef NORM
    #undef GRAD2GAP
    #undef FREE
    #undef CANCEL
}

/*  PG_TV2

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda ||x_i - x_(i-1)||_2 .
        
    To do so a Projected Gradient algorithm is used to solve its dual problem.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
*/
int PG_TV2(double *y,double lambda,double *x,double *info,int n){
    int nn=n-1,i;
    double stop,tmp,tmp2,lam,pNorm,qNorm,pNormSq,dist;
    double *p,*aux;
    lapack_int one=1,rc,nnp=nn;
    
    /* Macros */
            
    #define GRAD2GAP(w,g,gap,lambda,n,i,tmp) \
        gap = tmp = 0; \
        for(i=0;i<n;i++){ \
            tmp += g[i]*g[i]; \
            gap += w[i] * g[i]; \
        } \
        gap += lambda * sqrt(tmp); \
        gap = fabs(gap);
        
    #define NORM(x,n,i,tmp) \
        tmp = 0; \
        for(i=0;i<n;i++) tmp += x[i]*x[i]; tmp = sqrt(tmp);
    
    #define FREE \
        if(p) free(p); \
        if(aux) free(aux);
        
    #define CANCEL(txt,info) \
        printf("more_TV2: %s\n",txt); \
        FREE \
        info[INFO_RC] = RC_ERROR;\
        return 0;
        
    #define STEP 0.25
    #define MAX_PG 100000
    
    /* Alloc memory */
    p = (double*)calloc(nn,sizeof(double));
    aux = (double*)malloc(sizeof(double)*nn);
    if(!p || !aux){CANCEL("out of memory",info)}
    
    stop = DBL_MAX;
    info[INFO_ITERS] = 0;
    
    /* Construct problem */
    for(i=0;i<nn;i++)
        aux[i] = y[i] - y[i+1];
    
    /* Projected Gradient iterations */
    while(stop > STOP_MS && info[INFO_ITERS] < MAX_PG){
        /* Gradient step */
        for(i=0;i<nn;i++) p[i] = p[i] - STEP * aux[i];
        
        /* Projection step */
        NORM(p,n,i,tmp)
        if(tmp > lambda){
            tmp = lambda / tmp;
            for(i=0;i<nn;i++) p[i] *= tmp;
        }
        
        DUAL2PRIMAL(p,x,i)
        PRIMAL2GRAD(x,aux,i)
        GRAD2GAP(p,aux,stop,lambda,nn,i,tmp)

        info[INFO_ITERS]++;
    }

    info[INFO_GAP] = stop;
    
    /* Termination check */
    if(info[INFO_ITERS] >= MAX_PG){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(PG_TV2) WARNING: maximum number of iterations reached (%d).\n",MAX_PG);
        #endif
        info[INFO_RC] = RC_ITERS;
    }
    else info[INFO_RC] = RC_OK;
    
    FREE
    return 1;
    
    #undef NORM
    #undef FREE
    #undef CANCEL
    #undef STEP
    #undef MAX_PG
    #undef GRAD2GAP
}


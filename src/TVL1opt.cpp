/**
    Optimizers for problems dealing with TV-L1 norm regularization.
    
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

/*  PN_TV1

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda sum_i |x_i - x_(i-1)| .
        
    To do so a Projected Newton algorithm is used to solve its dual problem.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - info: array in which to store optimizer information.
        - n: length of array y (and x).
        - sigma: tolerance for sufficient descent.
        - ws: workspace of allocated memory to use. If NULL, any needed memory is locally managed.
*/
int PN_TV1(double *y,double lambda,double *x,double *info,int n,double sigma,Workspace *ws){
    int i,ind,nI,recomp,found,iters,nn=n-1;
    double lambdaMax,tmp,fval0,fval1,gRd,delta,grad0,stop,stopPrev,improve,rhs,maxStep,prevDelta;
    double *w=NULL,*g=NULL,*d=NULL,*aux=NULL,*aux2=NULL;
    int *inactive=NULL;
    lapack_int one=1,rc,nnp=nn,nIp;
    
    /* Macros */
            
    #define GRAD2GAP(g,w,gap,i) \
        gap = 0; \
        for(i=0;i<nn;i++) \
            gap += fabs(g[i]) * lambda + w[i] * g[i];
        
    #define PRIMAL2VAL(x,val,i) \
        val = 0; \
        for(i=0;i<n;i++) \
            val += x[i]*x[i]; \
        val *= 0.5;
            
    #define PROJECTION(w) \
        for(i=0;i<nn;i++) \
            if(w[i] > lambda) w[i] = lambda; \
            else if(w[i] < -lambda) w[i] = -lambda;
            
    #define CHECK_INACTIVE(w,g,inactive,nI,i) \
        for(i=nI=0 ; i<nn ; i++) \
            if( (w[i] > -lambda && w[i] < lambda) || (w[i] == -lambda && g[i] < -EPSILON) || (w[i] == lambda && g[i] > EPSILON) )  \
                inactive[nI++] = i;
                
    #define FREE \
        if(!ws){ \
            if(w) free(w); \
            if(g) free(g); \
            if(d) free(d); \
            if(aux) free(aux); \
            if(aux2) free(aux2); \
            if(inactive) free(inactive); \
        }
        
    #define CANCEL(txt,info) \
        printf("PN_TV1: %s\n",txt); \
        FREE \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;
            
    /* Alloc memory if no workspace available */
    if(!ws){
        w = (double*)malloc(sizeof(double)*nn);
        g = (double*)malloc(sizeof(double)*nn);
        d = (double*)malloc(sizeof(double)*nn);
        aux = (double*)malloc(sizeof(double)*nn);
        aux2 = (double*)malloc(sizeof(double)*nn);
        inactive = (int*)malloc(sizeof(int)*nn);
    }
    /* If a workspace is available, request memory */
    else{
        w = getDoubleWorkspace(ws);
        g = getDoubleWorkspace(ws);
        d = getDoubleWorkspace(ws);
        aux = getDoubleWorkspace(ws);
        aux2 = getDoubleWorkspace(ws);
        inactive = getIntWorkspace(ws);
    }
    if(!w || !g || ! d || !aux || !aux2 || !inactive)
        {CANCEL("out of memory",info)}

    /* Precompute useful quantities */
    for(i=0;i<nn;i++)
        w[i] = y[i+1] - y[i]; /* Dy */
    iters = 0;
        
    /* Factorize Hessian */
    for(i=0;i<nn-1;i++){
        aux[i] = 2;
        aux2[i] = -1;
    }
    aux[nn-1] = 2;
    dpttrf_(&nnp,aux,aux2,&rc);
    /* Solve Choleski-like linear system to obtain unconstrained solution */
    dpttrs_(&nnp, &one, aux, aux2, w, &nnp, &rc);
    
    /* Compute maximum effective penalty */
    lambdaMax = 0;
    for(i=0;i<nn;i++) 
        if((tmp = fabs(w[i])) > lambdaMax) lambdaMax = tmp;

    /* Check if the unconstrained solution is feasible for the given lambda */
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"lambda=%lf,lambdaMax=%lf\n",lambda,lambdaMax);
    #endif
    if(lambda >= lambdaMax){
        /* In this case all entries of the primal solution should be the same as the mean of y */
        tmp = 0;
        for(i=0;i<n;i++) tmp += y[i];
        tmp /= n;
        for(i=0;i<n;i++) x[i] = tmp;
        /* Gradient evaluation */
        PRIMAL2GRAD(x,g,i)
        /* Compute dual gap */
        GRAD2GAP(g,w,stop,i)
        if(info){
            info[INFO_GAP] = fabs(stop);
            info[INFO_ITERS] = 0;
            info[INFO_RC] = RC_OK;
        }
        FREE
        return 1;
    }
    
    /* If restart info available, use it to decide starting point */
    if(ws && ws->warm)
        memcpy(w,ws->warmDual,sizeof(double)*(n-1));
    
    /* Initial guess and gradient */
    PROJECTION(w)
    DUAL2PRIMAL(w,x,i)
    PRIMAL2GRAD(x,g,i)
    PRIMAL2VAL(x,fval0,i)

    /* Identify inactive constraints at the starting point */
    CHECK_INACTIVE(w,g,inactive,nI,i)
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"---------Starting point--------\n");
        fprintf(DEBUG_FILE,"w=["); for(i=0;i<nn;i++) fprintf(DEBUG_FILE,"%lf ",w[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"g=["); for(i=0;i<nn;i++) fprintf(DEBUG_FILE,"%lf ",g[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"inactive=["); for(i=0;i<nI;i++) fprintf(DEBUG_FILE,"%d ",inactive[i]); fprintf(DEBUG_FILE,"]\n");
        fprintf(DEBUG_FILE,"fVal=%lf\n",fval0);
        fprintf(DEBUG_FILE,"-------------------------------\n");
    #endif
    
    /* Solver loop */
    stop = DBL_MAX;
    stopPrev = -DBL_MAX;
    iters = 0;
    while(stop > STOP_PN && iters < MAX_ITERS_PN && fabs(stop-stopPrev) > EPSILON){
        /* If every constraint is active, we have finished */
        if(!nI){ FREE return 1;}
        
        /* Compute reduced Hessian (only inactive rows/columns) */
        for(i=0;i<nI-1;i++){
            aux[i] = 2;
            if(inactive[i+1]-inactive[i]!=1) aux2[i] = 0;
            else aux2[i] = -1;
        }
        aux[i] = 2;
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"alpha=["); for(i=0;i<nI;i++) fprintf(DEBUG_FILE,"%lf ",aux[i]); fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"beta=["); for(i=0;i<nI-1;i++) fprintf(DEBUG_FILE,"%lf ",aux2[i]); fprintf(DEBUG_FILE,"]\n");
        #endif
        /* Factorize reduced Hessian */
        nIp = nI;
        dpttrf_(&nIp,aux,aux2,&rc);
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"c=["); for(i=0;i<nI;i++) fprintf(DEBUG_FILE,"%lf ",aux[i]); fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"l=["); for(i=0;i<nI-1;i++) fprintf(DEBUG_FILE,"%lf ",aux2[i]); fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Solve Choleski-like linear system to obtain Newton updating direction */
        for(i=0;i<nI;i++)
            d[i] = g[inactive[i]];
        dpttrs_(&nIp, &one, aux, aux2, d, &nIp, &rc);
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"d=["); for(i=0;i<nI;i++) fprintf(DEBUG_FILE,"%lf ",d[i]); fprintf(DEBUG_FILE,"]\n");
        #endif
        
        /* Stepsize selection algorithm (quadratic interpolation) */
        gRd = 0;
        for(i=0;i<nI;i++)
            gRd += g[inactive[i]] * d[i];
        recomp = 0; delta = 1; found = 0;
        memcpy((void*)aux,(void*)w,sizeof(double)*nn);
        while(!found){
            /* Compute projected point after update */
            for(i=0;i<nI;i++){
                ind = inactive[i];
                aux[ind] = w[ind] - delta*d[i];
                if(aux[ind] > lambda) aux[ind] = lambda;
                else if(aux[ind] < -lambda) aux[ind] = -lambda;
            }
            /* Get primal point */
            DUAL2PRIMAL(aux,x,i)
            /* Compute new value of the objective function */
            PRIMAL2VAL(x,fval1,i)
            improve = fval0 - fval1;
            /* If zero improvement, the updating direction is not useful */
            if(improve <= EPSILON)
                break;
            /* Compute right hand side of Armijo rule */
            rhs = sigma * delta * gRd;
            /* Check if the rule is met */
            if(improve >= rhs) found = 1;
            else{
                if(!recomp){
                    /* Compute maximum useful stepsize */
                    maxStep = -DBL_MAX;
                    for(i=0;i<nI;i++){
                        if(d[i] < 0){
                            if((tmp=(w[inactive[i]]-lambda)/d[i]) > maxStep) maxStep = tmp;
                        } else if(d[i] > 0 && (tmp=(w[inactive[i]]+lambda)/d[i]) > maxStep) maxStep = tmp;
                    }
                    #ifdef DEBUG
                        fprintf(DEBUG_FILE,"maxStep=%lf\n",maxStep);
                    #endif
                    
                    /* Compute gradient w.r.t stepsize at the present position */
                    grad0 = 0;
                    if(!inactive[0]){
                        if( !((lambda == w[0]) && d[0] > 0) && !((lambda == w[0]) && d[0] < 0) )
                            grad0 += -d[0] * (2*w[0] - w[1] - y[1] + y[0]);
                    }
                    for(i=1;i<nI-1;i++){
                        ind = inactive[i];
                        if( !((lambda == w[ind]) && d[i] > 0) && !((lambda == w[ind]) && d[i] < 0) )
                            grad0 += -d[i] * (2*w[ind] - w[ind+1] - w[ind-1] - y[ind+1] + y[ind]);
                    }
                    if(inactive[nI-1] == nn-1){
                        if( !((lambda == w[nn-1]) && d[nI-1] > 0) && !((lambda == w[nn-1]) && d[nI-1] < 0) )
                            grad0 += -d[nI-1] * (2*w[nn-1] - w[nn-2] - y[nn] + y[nn-1]);
                    }
                
                    recomp = 1;
                }
                /* Use quadratic interpolation to determine next stepsize */
                tmp = grad0 * delta;
                prevDelta = delta;
                delta = - (tmp*delta) / (2 * (-improve - tmp));
                /* If larger than maximum stepsize, clip */
                if(delta > maxStep) delta = maxStep;
                /* If too similar to previous stepsize or larger, cut in half */
                if(delta-prevDelta >= -EPSILON) delta = prevDelta/2;
                /* If negative or zero, stop! */
                if(delta < EPSILON) found = true;
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"delta=%lf\n",delta);
                #endif
                /* Readjust maximum allowed step */
                maxStep = delta;
            }
        }
            
        /* Perform update */
        memcpy((void*)w,(void*)aux,sizeof(double)*nn);
        fval0 = fval1;
        
        /* Reconstruct gradient */
        PRIMAL2GRAD(x,g,i)

        /* Identify active and inactive constraints */
        CHECK_INACTIVE(w,g,inactive,nI,i)

        /* Compute stopping criterion */
        stopPrev = stop;
        GRAD2GAP(g,w,stop,i)
        
        iters++;
        
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"---------End of iteration %d--------\n",iters);
            fprintf(DEBUG_FILE,"w=["); for(i=0;i<nn;i++) fprintf(DEBUG_FILE,"%lf ",w[i]); fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"g=["); for(i=0;i<nn;i++) fprintf(DEBUG_FILE,"%lf ",g[i]); fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"inactive=["); for(i=0;i<nI;i++) fprintf(DEBUG_FILE,"%d ",inactive[i]); fprintf(DEBUG_FILE,"]\n");
            fprintf(DEBUG_FILE,"fVal=%lf\n",fval0);
            fprintf(DEBUG_FILE,"stop=%lf\n",stop);
        #endif
    }
    
    /* Termination check */
    if(iters >= MAX_ITERS_PN){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(PN_TV1) WARNING: maximum number of iterations reached (%d).\n",MAX_ITERS_PN);
        #endif
        if(info)
            info[INFO_RC] = RC_ITERS;
    }
    else if(fabs(stop-stopPrev) <= EPSILON){
        #ifdef DEBUG
            fprintf(DEBUG_FILE,"(PN_TV1) WARNING: search stuck, improvement is not possible.\n");
        #endif
        if(info)
            info[INFO_RC] = RC_STUCK;
    }
    else if(info) info[INFO_RC] = RC_OK;

    if(info){
        info[INFO_ITERS] = iters;
        info[INFO_GAP] = fabs(stop);
    }
    
    /* If restart structure available, store info for later warm restart */
    if(ws){
        memcpy(ws->warmDual,w,sizeof(double)*(n-1));
        ws->warm = 1;
    }
    
    FREE
    return 1;
         
    #undef GRAD2GAP        
    #undef PRIMAL2VAL            
    #undef PROJECTION            
    #undef CHECK_INACTIVE
    #undef FREE
    #undef CANCEL
}

/*  tautString_TV1

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator
    
        min_x 0.5 ||x-y||^2 + lambda sum_i |x_i - x_(i-1)| .
        
    To do so a Taut String algorithm is used to solve its equivalent problem.
    
    Inputs:
        - y: reference signal.
        - lambda: penalty parameter.
        - x: array in which to store the solution.
        - n: length of array y (and x).
*/
int tautString_TV1(double *y,double lambda,double *x,int n) {
    /* Minorant and minorant slopes */
    double mn, mx;
    /* Relative height of majorant and minorant slopes at the current points w.r.t. the tube center */
    double mnHeight, mxHeight;
    /* Last break points of minorant and majorant */
    int mnBreak, mxBreak;
    /* Last break point of taut string */
    int lastBreak;
    /* Auxiliary variables */
    int i, j;
    /* Helpful constants */
    const double minuslambda = -lambda;
    const double lambda2 = 2*lambda;
    const double minuslambda2 = 2*minuslambda;
        
    #define CANCEL(txt,info) \
        printf("tautString_TV1: %s\n",txt); \
        return 0;
    
    /* Starting point */
    mnHeight = mxHeight = 0;
    mn = minuslambda + y[0];
    mx = lambda + y[0];
    lastBreak = -1;
    mnBreak = mxBreak = 0;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Starting taut-string with length=%d and penalty=%lf\n",n,lambda); fflush(DEBUG_FILE);
    #endif
        
    /* Proceed along string */
    i = 0;
    while ( i < n ) {
        /* Loop over all points except the last one, that needs special care */
        while ( i < n-1 ) {
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mnHeight += mn - y[i];
        
            /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
            /* Majorant is r + lambda (except for last point), which is computed on the fly */   
            if ( lambda < mnHeight ) {
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,lambda,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
                #endif
                /* Break segment at last minorant breaking point */
                i = mnBreak + 1;
                /* Build valid segment up to this point using the minorant slope */
                for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
                    x[j] = mn;
                /* Start new segment after the break */
                lastBreak = mnBreak;
                /* Build first point of new segment, which can be done in closed form */
                mn = y[i]; 
                mx = lambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++;
                continue;
            }
            
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mxHeight += mx - y[i];
            
            /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
            /* Minorant is r - lambda (except for last point), which is computed on the fly */
            if ( minuslambda > mxHeight ) {
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,minuslambda,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
                #endif
                /* If violated, break segment at last majorant breaking point */
                i = mxBreak + 1;
                /* Build valid segment up to this point using the majorant slope */
                for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
                    x[j] = mx;
                /* Start new segment after the break*/
                lastBreak = mxBreak;
                /* Build first point of new segment, which can be done in closed form */
                mx = y[i]; 
                mn = minuslambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++;
                continue;
            }
            
            /* No violations at this point */

            /* Check if proyected majorant height is above ceiling */
            if ( mxHeight >= lambda ) {
                /* Update majorant slope */
                mx += ( lambda - mxHeight ) / ( i - lastBreak );
                /* Get correct majorant height (we are touching it!) */
                mxHeight = lambda;
                /* This is a possible majorant breaking point */
                mxBreak = i;
            }
            
            /* Check if proyected minorant height is under actual minorant */
            if ( mnHeight <= minuslambda ) {
                /* Update minorant slope */
                mn += ( minuslambda - mnHeight ) / ( i - lastBreak );
                /* Compute correct minorant height (we are touching it!) */
                mnHeight = minuslambda;
                /* This is a possible minorant breaking point */
                mnBreak = i;
            }

            #ifdef DEBUG
                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            
            /* At this point: no violations, so keep up building current segment */
            i++;
        }
        
        /* Special case i == n-1 (last point) */
        /* We try to validate the last segment, and if we can, we are finished */
        /* The code is essentially the same as the one for the general case, 
           the only different being that here the tube ceiling and floor are both 0 */
        
        /* Update height of minorant slope w.r.t. tube center */
        /* This takes into account both the slope of the minorant and the change in the tube center */
        mnHeight += mn - y[i];
    
        /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
        /* Majorant is 0 at this point */   
        if ( IS_POSITIVE(mnHeight) ) { // 0 < mnHeight
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"ENDING CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            /* Break segment at last minorant breaking point */
            i = mnBreak + 1;
            /* Build valid segment up to this point using the minorant slope */
            for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
                x[j] = mn;
            /* Start new segment after the break */
            lastBreak = mnBreak;
            /* Go back to main loop, starting a new segment */
            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
            mn = y[i]; 
            mx = lambda2+y[i];
            mxHeight = mnHeight = minuslambda;
            mnBreak = mxBreak = i;
            continue;
        }
            
        /* Update height of minorant slope w.r.t. tube center */
        /* This takes into account both the slope of the minorant and the change in the tube center */
        mxHeight += mx - y[i];
        
        /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
        /* Minorant is 0 at this point */
        if ( IS_NEGATIVE(mxHeight) ) { // 0 > mxHeight
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"ENDING FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            /* If violated, break segment at last majorant breaking point */
            i = mxBreak + 1;
            /* Build valid segment up to this point using the majorant slope */
            for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
                x[j] = mx;
            /* Start new segment after the break*/
            lastBreak = mxBreak;
            /* Go back to main loop, starting a new segment */
            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
            mx = y[i]; 
            mn = minuslambda2+y[i];
            mxHeight = mnHeight = lambda;
            mnBreak = mxBreak = i;
            continue;
        }
        
        /* No violations at this point */
        
        /* Check if proyected minorant height is under actual minorant */
        if ( mnHeight <= 0 ) {
            /* Update minorant slope */
            mn += ( - mnHeight ) / ( i - lastBreak );
        }

        #ifdef DEBUG
            fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
        #endif
        
        /* At this point: we are finished validating last segment! */
        i++;
    }
    
    /* Build last valid segment */
    for ( i = lastBreak+1 ; i < n ; i++ )
        x[i] = mn;
        
    /* Return */
    return 1;
    
    #undef CANCEL
}


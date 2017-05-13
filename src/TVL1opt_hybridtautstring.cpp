/**
    Hybrid taut-string solver for TV-L1 proximity
    
    @author Álvaro Barbero Jiménez
*/
#include <math.h>
#include <stdio.h>
#include "TVopt.h"

// Default exponent for maximum backtracking steps
#define BACKTRACKSEXP 1.05

/*
    Checks the current number of backtracks, and if sufficently large
    switches to concave-convex taut-string method.

    Arguments:
        - backtracks: number of backtracked values in the algorithm so far
        - maxbacktracks: maximum number of backtrack steps to perform
        - signal: pointer to input signal
        - n: length of input signal
        - lam: regularization parameter
        - prox: pointer to solution built so far
        - currentstep: value in the input up to a solution has been found
        - offset: y offset of the current point with respect to the tube center

    The macro has no effect if no switch was performed. If it was,
    forces immediate return of the method, and the prox array contains the full
    solution to the problem. 
*/
#define ATTEMPTSWITCH(backtracks, maxbacktracks, signal, n, lam, prox, currentstep, offset) \
    if (backtracks > maxbacktracks && currentstep + 1 < n) { \
        classicTautString_TV1_offset(signal+currentstep, n-currentstep, lam, prox+currentstep, offset); \
        return; \
    }

/*
    Solves the Total-Variation l1 problem through a hybrid method
    that selects the best approach for each case.

    The method starts by running Condat's method (linearized taut-string)
    and if the amount of backtracking starts being too high we
    switch to the classic (concave-convex) taut-string method.

    Arguments:
        - y: input signal
        - n: length of signal
        - lambda: strength of l1 regularization
        - x: array in which to store solution
        - backtracksexp: maximum number of backtrack steps, defined as an 
        exponent over the number of elements in the signal. Maximum number 
        is 2 as the worst case for Condat's method is O(n^2)
        
    Returns: signal after TV-l1 proximity.
*/
void hybridTautString_TV1_custom(double *y, int n, double lambda, double *x, double backtracksexp) {
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
    // Backtracking cost tracker
    int backtracks = 0;
    double maxbacktracks = pow(n,backtracksexp);
    
    /* Starting point */
    mnHeight = mxHeight = 0;
    mn = minuslambda + y[0];
    mx = lambda + y[0];
    lastBreak = -1;
    mnBreak = mxBreak = 0;
        
    /* Proceed along string */
    i = 0;
    while ( i < n ) {
        /* Loop over all points except the last one, that needs special care */
        while ( i < n-1 ) {
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mnHeight += mn - y[i];
        
            /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
            /* Majorant is r + lambda (except for last point), which is computed on the fly */   
            if ( lambda < mnHeight ) {
                /* Break segment at last minorant breaking point */
                i = mnBreak + 1;
                /* Build valid segment up to this point using the minorant slope */
                for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
                    x[j] = mn;
                /* Test if we must switch to convex-concave method */
                ATTEMPTSWITCH(backtracks, maxbacktracks, y, n, lambda, x, i, -lambda)
                /* Start new segment after the break */
                lastBreak = mnBreak;
                /* Build first point of new segment, which can be done in closed form */
                mn = y[i]; 
                mx = lambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++; backtracks++;
                continue;
            }
            
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mxHeight += mx - y[i];
            
            /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
            /* Minorant is r - lambda (except for last point), which is computed on the fly */
            if ( minuslambda > mxHeight ) {
                /* If violated, break segment at last majorant breaking point */
                i = mxBreak + 1;
                /* Build valid segment up to this point using the majorant slope */
                for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
                    x[j] = mx;
                /* Test if we must switch to convex-concave method */
                ATTEMPTSWITCH(backtracks, maxbacktracks, y, n, lambda, x, i, lambda)
                /* Start new segment after the break*/
                lastBreak = mxBreak;
                /* Build first point of new segment, which can be done in closed form */
                mx = y[i]; 
                mn = minuslambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++; backtracks++;
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
            
            /* At this point: no violations, so keep up building current segment */
            i++; backtracks++;
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
            /* Break segment at last minorant breaking point */
            i = mnBreak + 1;
            /* Build valid segment up to this point using the minorant slope */
            for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
                x[j] = mn;
            /* Test if we must switch to convex-concave method */
            ATTEMPTSWITCH(backtracks, maxbacktracks, y, n, lambda, x, i, -lambda)
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
            /* If violated, break segment at last majorant breaking point */
            i = mxBreak + 1;
            /* Build valid segment up to this point using the majorant slope */
            for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
                x[j] = mx;
            /* Test if we must switch to convex-concave method */
            ATTEMPTSWITCH(backtracks, maxbacktracks, y, n, lambda, x, i, lambda)
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
        
        /* At this point: we are finished validating last segment! */
        i++;
    }
    
    /* Build last valid segment */
    for ( i = lastBreak+1 ; i < n ; i++ )
        x[i] = mn;
}

void hybridTautString_TV1(double *y, int n, double lambda, double *x) {
    hybridTautString_TV1_custom(y, n, lambda, x, BACKTRACKSEXP);
}

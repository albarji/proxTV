/**
    General functions dealing with TV norm regularization.
    
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

/**
    General function for 1 dimensional unweighted TV optimization. Invokes different algorithms depending on the
    TV norm p.
    
    @param y input signal
    @param lambda regularization value
    @param x preallocated space in which to store solution
    @param info preallocated space in which to store optimizer information
    @param n length of input signal
    @param p degree of the TV norm to apply
    @param ws pointer to workspace to use, or null if the method should alloc its own memory
*/
int TV(double *y,double lambda,double *x,double *info,int n,double p,Workspace *ws) {
    #define CANCEL(txt,info) \
        printf("TVopt: %s\n",txt); \
        if(info) info[INFO_RC] = RC_ERROR;\
        return 0;

    // Check arguments
    if ( p < 1 )
        {CANCEL("TV only works for norms p >= 1",info)}
        
    // Invoke appropriate method
    if(p == 1) {
        tautString_TV1(y, lambda, x, n);
        if ( info ) {
            info[INFO_RC] = RC_OK;
            info[INFO_ITERS] = 0;
            info[INFO_GAP] = 0;
        }
    }
    else if(p == 2)
        morePG_TV2(y, lambda, x, info, n, ws);
    else // general p norm
        GPFW_TVp(y, lambda, x, info, n, p, ws); 
        
    return 1;
 
    #undef CANCEL   
}

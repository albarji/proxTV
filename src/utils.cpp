/**
    General utility functions.
    
    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "utils.h"

/**
    Returns the sign of a number.
    
    @argument s input number
    @returns -1 if number is negative, +1 if positive, 0 if exactly 0
*/
short sign(double s) {
    return ( s >= 0 ) - ( s <= 0 );
}

/**
    Returns the minimum of two numbers
    
    @argument x
    @argument y
    @returns value of the minimum number
*/
double min(double x, double y) {
    if ( x < y )
        return x;
    else 
        return y;
}

/**
    Returns the maximum of two numbers
    
    @argument x
    @argument y
    @returns value of the maximum number
*/
double max(double x, double y) {
    if ( x > y )
        return x;
    else 
        return y;
}

/**
    Performs radial projection of a given point x.
    The present norm of x and desired norm lambda must be provided.

    @argument x point to project.
    @argument n length of x.
    @argument norm present norm of x.
    @argument lambda desired norm of x.
*/
void radialProjection(double *x, int n, double norm, double lambda) {
    /* Check projection is needed */
    if ( norm > lambda ) {
        /* Compute scaling factor */
        double factor = lambda / norm;
        /* Scale all components */
        for (int i = 0 ; i < n ; i++ )
            x[i] *= factor;
    }
}

/**
    Creates and initiallizes a workspace structure
    
    @argument n dimension of data vectors.
    @returns pointer to new workspace
*/
Workspace* newWorkspace(int n) {
    Workspace *ws;
    int i;
    
    #define CANCEL \
        freeWorkspace(ws); \
        return NULL;

    /* Alloc structure */
    ws = (Workspace*)calloc(1,sizeof(Workspace));
    if(!ws) {CANCEL}
    
    /* Alloc input and output fields */
    ws->n = n;
    ws->in = (double*)malloc(sizeof(double)*n);
    ws->out = (double*)malloc(sizeof(double)*n);
    if(!ws->in || !ws->out) {CANCEL}
    
    /* Alloc generic double fields holder */
    ws->d = (double**)calloc(WS_MAX_MEMORIES,sizeof(double*));
    if(!ws->d) {CANCEL}
    ws->nd = ws->maxd = 0;
    
    /* Alloc generic int fields holder */
    ws->i = (int**)calloc(WS_MAX_MEMORIES,sizeof(int*));
    if(!ws->i) {CANCEL}
    ws->ni = ws->maxi = 0;
    
    /* Alloc warm restart fields */
    ws->warmDual = (double*)malloc(sizeof(double)*(n-1));
    if(!ws->warmDual) {CANCEL}
        
    return ws;
    
    #undef CANCEL
}

/**
    Resets the memory assignments in the workspace
    
    @argument ws pointer to workspace to reset
*/
void resetWorkspace(Workspace *ws) {
    /* Reset assigned memory counters */
    ws->ni = ws->nd = 0;
}

/**
    Returns a vector of double memory from the workspace.
    New memory is allocated on-demand
    
    @argument ws pointer to workspace
*/
double* getDoubleWorkspace(Workspace *ws) {
    /* Check if maximum memory limit has been topped */
    if (ws->nd == WS_MAX_MEMORIES)
        return NULL;

    /* Check if additional memory allocation is necessary */
    if (ws->nd == ws->maxd) {
        /* Alloc additional memory */
        ws->d[ws->nd] = (double*)calloc(ws->n,sizeof(double));
        if (!ws->d[ws->nd])
            return NULL;
        /* Increase number of allocated vectors */
        ws->maxd++;
    }
    
    /* Increase number of used vectors */
    ws->nd++;
    
    /* Return memory */
    return ws->d[ws->nd-1];
}

/**
    Returns a vector of int memory from the workspace.
    New memory is allocated on-demand
    
    @argument ws pointer to workspace
*/
int* getIntWorkspace(Workspace *ws) {
    /* Check if maximum memory limit has been topped */
    if (ws->ni == WS_MAX_MEMORIES)
        return NULL;

    /* Check if additional memory allocation is necessary */
    if (ws->ni == ws->maxi) {
        /* Alloc additional memory */
        ws->i[ws->ni] = (int*)calloc(ws->n,sizeof(int));
        if (!ws->i[ws->ni])
            return NULL;
        /* Increase number of allocated vectors */
        ws->maxi++;
    }
    
    /* Increase number of used vectors */
    ws->ni++;
    
    /* Return memory */
    return ws->i[ws->ni-1];
}

/* Frees a workspace structure */
void freeWorkspace(Workspace *ws){
    int i;

    if(ws){
        /* Input/output fields */
        if(ws->in) free(ws->in);
        if(ws->out) free(ws->out);
        /* Generic memory fields */
        if(ws->d){
            for(i=0;i<ws->nd;i++) if(ws->d[i]) free(ws->d[i]);
            free(ws->d);
        }
        if(ws->i){
            for(i=0;i<ws->ni;i++) if(ws->i[i]) free(ws->i[i]);
            free(ws->i);
        }
        /* Warm restart fields */
        if(ws->warmDual) free(ws->warmDual);
        
        free(ws);
    }
}

/* Allocs memory for an array of p workspaces */
Workspace** newWorkspaces(int n,int p){
    int i;
    Workspace **wa=NULL;
    
    #define CANCEL \
        freeWorkspaces(wa,p); \
        return NULL;
    
    /* Alloc the array */
    wa = (Workspace**)calloc(p,sizeof(Workspace*));
    if(!wa) {CANCEL}
    
    /* Alloc each individual workspace */
    for(i=0;i<p;i++){
        wa[i] = newWorkspace(n);
        if(!wa[i]) {CANCEL}
    }
    
    return wa;
}

/* Frees an array of p workspaces */
void freeWorkspaces(Workspace **wa,int p){
    int i;
    
    if(wa){
        /* Free each workspace */
        for(i=0;i<p;i++) freeWorkspace(wa[i]);
        /* Free array */
        free(wa);
    }
}

/* Comparison of floating point numbers, in reversed sense (to allow for descending sorting) */
int compareDoublesReversed(const void *v1, const void *v2){
    static double d1,d2;
    d1 = *((double*)v1);
    d2 = *((double*)v2);
    if(d1 < d2) return 1;
    else if(d1 > d2) return -1;
    else return 0;
}

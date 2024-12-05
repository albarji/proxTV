/**
    General utility functions.

    @author Álvaro Barbero Jiménez
    @author Suvrit Sra
*/

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/*** Math functions ***/
short sign(double s) {
    return ( s >= 0 ) - ( s <= 0 );
}

double min(double x, double y) {
    return (x < y) ? x : y;
}

double max(double x, double y) {
    return (x > y) ? x : y;
}

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

/*** Memory allocation functions ***/
Workspace* newWorkspace(int n) {
    Workspace *ws;

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

void resetWorkspace(Workspace *ws) {
    /* Reset assigned memory counters */
    ws->ni = ws->nd = 0;
}

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

void freeWorkspaces(Workspace **wa,int p){
    int i;

    if(wa){
        /* Free each workspace */
        for(i=0;i<p;i++) freeWorkspace(wa[i]);
        /* Free array */
        free(wa);
    }
}

int compareDoublesReversed(const void *v1, const void *v2){
    static double d1,d2;
    d1 = *((double*)v1);
    d2 = *((double*)v2);
    if(d1 < d2) return 1;
    else if(d1 > d2) return -1;
    else return 0;
}

#ifdef __cplusplus
}
#endif

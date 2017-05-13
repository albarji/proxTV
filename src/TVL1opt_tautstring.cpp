#include <stdlib.h>
#include <string.h>
#include "TVopt.h"

/********************************************
    Data Structures
********************************************/

/*
    Structure defining a bidimensional point with integer x value
*/
typedef struct {
    int x;
    double y;
} Point;

/*
    Structure defining a segment in a polyline

    The segment is represented by its span in the (x,y) dimension
    and its slope. This is an efficient representation for calculations
    in the taut string method, even if one of the quantities is redundant.
*/
typedef struct {
    int incx;
    double incy;
    double slope;
} Segment;

/*
    Structure defining a linear buffer of segments.

    The structure is an array of Segment, allowing insertion at
    the end and removal at both ends. It also allows constant
    time random access.

    To achieve this the structure must an amount of memory as
    large as the number of Points to be introduced.
*/
typedef struct {
    Segment* segments; // Array of segments in the buffer
    Segment* first; // Pointer to first valid segment in buffer
    Segment* last; // Pointer to last valid segment in buffer
} Buffer;

/*
    Creates a new Buffer of the given size
*/
inline Buffer* newPB(int n) {
    // Create structure, init as empty
    Buffer* pf = (Buffer*)calloc(1, sizeof(Buffer));
    // Create array of Segments
    pf->first = pf->segments = (Segment*)malloc(sizeof(Segment)*n);
    pf->last = pf->first-1;
    return pf;
}

/*
    Frees all memory reserved by a Buffer
*/
inline void freePB(Buffer* pb) {
    // Free array of segments
    free(pb->segments);
    // Free structure
    free(pb);
}

/*
    Adds a new Segment to a Buffer

    The values of the Segment are copied

    Arguments:
        Buffer* pb: pointer to Buffer from which to remove
        Segment* s: pointer to segment to add
*/
#define pushPB(pb, s) \
    *(++pb->last) = *(s);

/*
    Removes the first element in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer from which to remove
*/
#define removefirstPB(pb) \
    pb->first++;

/*
    Returns the first element in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer from which to read
*/
#define peekfirstPB(pb) \
    (pb->first)

/*
    Removes and returns the last element in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer from which to remove
        Segment* s: pointer to use to return last element
*/
#define poplastPB(pb, s) \
    s = pb->last--;

/*
    Returns the last element in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer from which to read
*/
#define peeklastPB(pb) \
    (pb->last)

/*
    Returns the number of Points in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer to measure
*/
#define sizePB(pb) \
    (pb->last - pb->first + 1)

/*
    Removes all elements in a Buffer

    Arguments:
        Buffer* pb: pointer to Buffer to clear
*/
#define clearPB(pb) \
    pb->first = pb->segments; \
    pb->last = pb->first-1;

/********************************************
    Taut-string functions
********************************************/

/*
    Adds a new segment to a concave majorant function

    Arguments:
        - majorant: Buffer* with majorant to augment
        - s: Segment* to add to majorant
        - saux: Segment* for auxiliary operations
        - iaux: int for auxiliary operations
*/
#define concavemajorantadd(majorant, s, saux, iaux) \
    if ( s->slope > peeklastPB(majorant)->slope )  { \
        iaux = sizePB(majorant); \
        do { \
            poplastPB(majorant, saux); \
            s->incx += saux->incx; \
            s->incy += saux->incy; \
            iaux--; \
        } while ( iaux >= 1 && s->incy > s->incx * peeklastPB(majorant)->slope ); \
        s->slope = s->incy / s->incx; \
    } \
    pushPB(majorant, s);

/*
    Adds a new point to a convex minorant function

    Arguments:
        - minorant: Buffer* with minorant to augment
        - s: Segment* to add to minorant
        - saux: Segment* for auxiliary operations
        - iaux: int for auxiliary operations
*/
#define convexminorantadd(minorant, s, saux, iaux) \
    if ( s->slope < peeklastPB(minorant)->slope )  { \
        iaux = sizePB(minorant); \
        do { \
            poplastPB(minorant, saux); \
            s->incx += saux->incx; \
            s->incy += saux->incy; \
            iaux--; \
        } while ( iaux >= 1 && s->incy < s->incx * peeklastPB(minorant)->slope ); \
        s->slope = s->incy / s->incx; \
    } \
    pushPB(minorant, s);

/*
    Returns a new taut-string knot for given crossing majorant/minorant
*/
inline Segment* newknot(Buffer* majorant, Buffer* minorant, Point* origin, Point* lastexplored, double lam) {
    // Get left-most segment in both approximators
    Segment *smin = peekfirstPB(minorant);
    Segment *smaj = peekfirstPB(majorant);
    Segment s;
    Segment *spointer = &s;
    // Shortest segment defines the new knot
    if ( smin->incx < smaj->incx ) {
        // Remove first segment, rest of the minorant is still valid
        removefirstPB(minorant);
        // Majorant is a single segment from new knot to last original point
        s.incx = lastexplored->x - origin->x - smin->incx;
        s.incy = lastexplored->y - lam - origin->y - smin->incy;
        s.slope = s.incy / s.incx;
        clearPB(majorant);
        pushPB(majorant, spointer);
        // Update origin
        origin->x += smin->incx;
        origin->y += smin->incy;
        return smin;
    }
    // Left-most majorant touching point is the new knot
    else {
        // Remove first point, rest of the majorant is still valid
        removefirstPB(majorant);
        // Minorant is a single segment from new knot to last original point
        s.incx = lastexplored->x - origin->x - smaj->incx;
        s.incy = lastexplored->y + lam - origin->y - smaj->incy;
        s.slope = s.incy / s.incx;
        clearPB(minorant);
        pushPB(minorant, spointer);
        // Update origin
        origin->x += smaj->incx;
        origin->y += smaj->incy;
        return smaj;
    }
}

/*
    Adds a new segment to the solution of the prox-TV problem

    Arguments:
        - prox: double* array in which to store prox-TV solution
        - segment: Segment* of the last fixed taut-string knot
        - j: int for auxiliary operations

    The new prox-TV segment is added to the prox array
*/
#define addproxsegment(prox, segment, j) \
    for ( j = 0 ; j < segment->incx ; j++ ) \
        (prox)[j] = segment->slope;

/*
    Solves the Total-Variation l1 problem

    To do so an efficient implementation of the classic taut-string method is used.

    Arguments:
        - signal: input signal
        - n: length of signal
        - lam: strength of l1 regularization
        - prox: array in which to save results
        - initial y offset in the start of the algorithm. Useful when
        switching to this algorithm from othe taut-string method.
        A positive offset mean the starting point is above the tube center
        by this offset quantity, negative means it is under the tube center.
        
    Returns: signal after TV-l1 proximity.
*/
int classicTautString_TV1_offset(double *signal, int n, double lam, double *prox, double offset) {
    // Deal with extreme cases (lam=0, n=1, etc...)
    if ( n <= 0 )
        return 1;
    if (lam <= 0 || n == 1 ) {
        memcpy(prox, signal, n*sizeof(double));
        return 1;
    }

    // Required memory structures
    Buffer* majorant = newPB(n);
    Buffer* minorant = newPB(n);

    // Start by adding the initial segment from (0,0) to the first point
    // Crossing checks are not necessary as crosses are impossible
    Segment segment;
    Segment *dirsegment = &segment;
    segment.incx = 1;
    segment.slope = segment.incy = signal[0] - offset - lam;
    pushPB(majorant, dirsegment);
    segment.slope = segment.incy = signal[0] - offset + lam;
    pushPB(minorant, dirsegment);

    // Initial point of current taut-string segment
    Point origin;
    Point *porigin = &origin;
    origin.x = 0;
    origin.y = offset;
    // Middle point of last explored point in the tube
    Point lastexplored;
    Point *plastexplored = &lastexplored;
    lastexplored.x = 1;
    lastexplored.y = signal[0];
        
    // Iterate along the signal length
    Segment *saux;
    int i, iaux;
    double *pwriter = prox;
    for ( i = 1 ; i < n-1 ; i++ ) {
        // Update majorant
        segment.incx = 1;
        segment.slope = segment.incy = signal[i];
        concavemajorantadd(majorant, dirsegment, saux, iaux);
        // Update minorant
        segment.incx = 1;
        segment.slope = segment.incy = signal[i];
        convexminorantadd(minorant, dirsegment, saux, iaux);
        // Update last explored point
        lastexplored.x++;
        lastexplored.y += signal[i];

        // Check for slope crossings at the first point
        while ( peekfirstPB(minorant)->slope < peekfirstPB(majorant)->slope ) {
            // Crossing detected!
            saux = newknot(majorant, minorant, porigin, plastexplored, lam);
            addproxsegment(pwriter, saux, iaux);
            pwriter += saux->incx;
        }
    }

    // Update majorant with last segment
    segment.incx = 1;
    segment.slope = segment.incy = signal[n-1] + lam;
    concavemajorantadd(majorant, dirsegment, saux, iaux);
    // Update minorant with last segment
    segment.incx = 1;
    segment.slope = segment.incy = signal[n-1] - lam;
    convexminorantadd(minorant, dirsegment, saux, iaux);

    // At this point, because the endpoint of the tube is the same
    // for both majorant and minorant, either the majorant or the minorant
    // is a single straight segment while the other can be larger.
    // The remaining of the taut-string must be the multi-segment component.
    Buffer *larger = (sizePB(majorant) > sizePB(minorant)) ? majorant : minorant;
    saux = peekfirstPB(larger);
    for ( i = sizePB(larger) ;  i > 0 ; i-- ) {
        addproxsegment(pwriter, saux, iaux);
        pwriter += (saux++)->incx;
    }

    // Free memory
    freePB(majorant);
    freePB(minorant);
    return 1;
}

/*
    Solves the Total-Variation l1 problem.

    To do so an efficient implementation of the classic taut-string method is used.

    Arguments:
        - signal: input signal
        - n: length of signal
        - lam: strength of l1 regularization
        - prox: array in which to save results
        
    Returns: signal after TV-l1 proximity.
*/
int classicTautString_TV1(double *signal, int n, double lam, double *prox) {
    return classicTautString_TV1_offset(signal, n, lam, prox, 0);
}


/**
    Lapack functions wrapping
    This file creates a NAME_wrap layer on top of lapack functions with name: NAME
    The wrap function will call lapack or eigen functions depending on the compilation options
    WITH LAPACK:
      Declare extern NAME_ functions that will be provided by LAPACKE
    WITH EIGEN:
      Declare NAME_eigen, equivalent to lapack but using EIGEN.

    @author Pablo Hernandez-Cerdan
*/
#ifdef NOMATLAB
#undef lapack_int
#define lapack_int int
#else // WITH_MATLAB
#include "mex.h"
#ifdef PROXTV_USE_LAPACK
#include "lapack.h"
#endif
#define lapack_int  ptrdiff_t
#endif

#ifdef PROXTV_USE_LAPACK
extern "C" {
    /**
     *  DPTTRS solves a tridiagonal system of the form
     *  A * X = B
     *  using the L*D*L' factorization of A computed by DPTTRF.  D is a
     *  diagonal matrix specified in the vector D, L is a unit bidiagonal
     *  matrix whose subdiagonal is specified in the vector E, and X and B
     *  are N by NRHS matrices.
     *
     *  Arguments
     *  =========
     *
     *  @param n (input) INTEGER
     *           The order of the tridiagonal matrix A.  N >= 0.
     *
     *  @param nrhs (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  @param d  (input) DOUBLE PRECISION array, dimension (N)
     *          The n diagonal elements of the diagonal matrix D from the
     *          L*D*L' factorization of A.
     *
     *  @param e  (input) DOUBLE PRECISION array, dimension (N-1)
     *          The (n-1) subdiagonal elements of the unit bidiagonal factor
     *          L from the L*D*L' factorization of A.  E can also be regarded
     *          as the superdiagonal of the unit bidiagonal factor U from the
     *          factorization A = U'*D*U.
     *
     *  @param b  (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side vectors B for the system of
     *          linear equations.
     *          On exit, the solution vectors, X.
     *
     *  @param ldb  (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  @param info (output) INTEGER
     *          = 0: successful exit
     *          < 0: if INFO = -k, the k-th argument had an illegal value
     *
     */
    void dpttrs_(lapack_int* n, lapack_int* nrhs, const double* d,
            const double* e, double* b, lapack_int* ldb, lapack_int *info );
    /**
     *  DPTTRF computes the L*D*L' factorization of a real symmetric
     *  positive definite tridiagonal matrix A.  The factorization may also
     *  be regarded as having the form A = U'*D*U.
     *
     *  Arguments
     *  =========
     *
     *  @param n (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  @param d (input/output) DOUBLE PRECISION array, dimension (N)
     *          On entry, the n diagonal elements of the tridiagonal matrix
     *          A.  On exit, the n diagonal elements of the diagonal matrix
     *          D from the L*D*L' factorization of A.
     *
     *  @param e (input/output) DOUBLE PRECISION array, dimension (N-1)
     *          On entry, the (n-1) subdiagonal elements of the tridiagonal
     *          matrix A.  On exit, the (n-1) subdiagonal elements of the
     *          unit bidiagonal factor L from the L*D*L' factorization of A.
     *          E can also be regarded as the superdiagonal of the unit
     *
     *  @param b  (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
     *          On entry, the right hand side vectors B for the system of
     *          linear equations.
     *          On exit, the solution vectors, X.
     *          bidiagonal factor U from the U'*D*U factorization of A.
     *
     *  @param info (output) INTEGER
     *          = 0: successful exit
     *          < 0: if INFO = -k, the k-th argument had an illegal value
     *          > 0: if INFO = k, the leading minor of order k is not
     *               positive definite; if k < N, the factorization could not
     *               be completed, while if k = N, the factorization was
     *               completed, but D(N) <= 0.
     */
    void dpttrf_( lapack_int* n, double* d, double* e, lapack_int *info );
}

#else // USE_EIGEN

/**
 *  DPTTRS solves a tridiagonal system of the form
 *  A * X = B
 *  using the L*D*L' factorization of A computed by DPTTRF.  D is a
 *  diagonal matrix specified in the vector D, L is a unit bidiagonal
 *  matrix whose subdiagonal is specified in the vector E, and X and B
 *  are N by NRHS matrices.
 *
 *  DPTTRF computes the L*D*L' factorization of a real symmetric
 *  positive definite tridiagonal matrix A.  The factorization may also
 *  be regarded as having the form A = U'*D*U.
 *
 *  Arguments
 *  =========
 *
 *  @param n (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  @param d (input) DOUBLE PRECISION array, dimension (N)
 *          On entry, the n diagonal elements of the tridiagonal matrix
 *          A. On exit, the n diagonal elements of the diagonal matrix
 *          D from the L*D*L' factorization of A.
 *
 *  @param e (input) DOUBLE PRECISION array, dimension (N-1)
 *          On entry, the (n-1) subdiagonal elements of the tridiagonal
 *          matrix A.  On exit, the (n-1) subdiagonal elements of the
 *          unit bidiagonal factor L from the L*D*L' factorization of A.
 *          E can also be regarded as the superdiagonal of the unit
 *
 *  @param b  (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the right hand side vectors B for the system of
 *          linear equations.
 *          On exit, the solution vectors, X.
 *
 *  @param info (output) INTEGER ----- UNUSED -----
 *          = 0: successful exit
 *          < 0: if INFO = -k, the k-th argument had an illegal value
 *          > 0: if INFO = k, the leading minor of order k is not
 *               positive definite; if k < N, the factorization could not
 *               be completed, while if k = N, the factorization was
 *               completed, but D(N) <= 0.
 */
void dpttrf_plus_dpttrs_eigen( lapack_int* n, double* d, double* e, double* b);

#endif

import os.path
from sys import platform as _platform

from cffi import FFI


ffi = FFI()
ffi.cdef("""
    typedef struct {
        ...;
    } Workspace;

    // Condat's implementation.
    void TV1D_denoise(double* input, double* output, const int width,
                      const double lambda);
    void TV1D_denoise_tautstring(double* input, double* output, int width, 
                                 const double lambda);
    // Ryan's implementation of Johnson's algorithm
    void dp(int n, double *y, double lam, double *beta);

    /* TV-L1 solvers */

    int PN_TV1(double *y,double lambda,double *x,double *info,int n,double sigma,Workspace *ws);
    int linearizedTautString_TV1(double *y,double lambda,double *x,int n);
    int classicTautString_TV1(double *signal, int n, double lam, double *prox);
    void hybridTautString_TV1(double *y, int n, double lambda, double *x);
    void hybridTautString_TV1_custom(double *y, int n, double lambda, double *x, double backtracksexp);
    void SolveTVConvexQuadratic_a1_nw(int n, double* b, double w, double* solution);

    /* Weighted TV-L1 solvers */
    int PN_TV1_Weighted(double* Y, double* W, double* X, double* info, int n,
                        double sigma, Workspace* ws);
    int tautString_TV1_Weighted(double *y, double* lambda, double *x, int n);
    void SolveTVConvexQuadratic_a1(int n, double* b, double* w, double* solution);

    /* TV-L2 solvers */
    int more_TV2(double *y,double lambda, double *x, double *info, int n);
    int PG_TV2(double *y, double lambda, double *x,double *info, int n);
    int morePG_TV2(double *y, double lambda, double *x, double *info, int n,
                   Workspace *ws);

    /* Weighted TV-L2 solvers */
    int DR2L1W_TV(size_t M, size_t N, double* unary, double*W1, double*W2,
                  double *s, int nThreads, int maxit, double* info);


    /* 2-dimensional TV solvers */
    int PD2_TV(double *y, double *lambdas, double *norms, double *dims,
               double *x, double *info, int *ns, int nds, int npen, int ncores,
               int maxIters);
    int DR2_TV(size_t M, size_t N, double*unary, double W1, double W2,
               double norm1, double norm2, double*s, int nThreads, int maxit,
               double* info);
    int CondatChambollePock2_TV(size_t M, size_t N, double*Y, double lambda,
                                double*X, short alg, int maxit, double* info);
    int Yang2_TV(size_t M, size_t N, double*Y, double lambda, double*X,
                 int maxit, double* info);
    int Kolmogorov2_TV(size_t M, size_t N, double*Y, double lambda, double*X, 
                 int maxit, double* info);

    /* General-dimension TV solvers */
    int PD_TV(double *y, double *lambdas, double *norms, double *dims,
                double *x, double *info, int *ns, int nds, int npen,
                int ncores, int maxIters);

    /* TV-Lp solvers */
    int GP_TVp(double *y, double lambda, double *x, double *info, int n,
               double p, Workspace *ws);
    int OGP_TVp(double *y, double lambda, double *x, double *info, int n,
                double p, Workspace *ws);
    int FISTA_TVp(double *y, double lambda, double *x, double *info, int n,
                  double p, Workspace *ws);
    int FW_TVp(double *y, double lambda, double *x, double *info, int n,
               double p, Workspace *ws);
    int GPFW_TVp(double *y, double lambda, double *x, double *info, int n,
                 double p, Workspace *ws);
""")

sources = [os.path.join('src', fname) for fname in (
    'condat_fast_tv.cpp', 'johnsonRyanTV.cpp', 'TVL1opt_kolmogorov.cpp', 
    'LPopt.cpp', 'TV2Dopt.cpp', 'TV2DWopt.cpp', 'TVgenopt.cpp', 'TVL1opt.cpp', 
    'TVL1opt_tautstring.cpp', 'TVL1opt_hybridtautstring.cpp', 'TVL1Wopt.cpp', 
    'TVL2opt.cpp', 'TVLPopt.cpp', 'TVNDopt.cpp', 'utils.cpp'
)]

extra_compile_args = []
extra_link_args = []
if _platform == 'darwin':
    # if openblas was installed by homerew is present use this for lapacke.h
    if os.path.exists('/usr/local/opt/openblas/include'):
        extra_compile_args.append('-I/usr/local/opt/openblas/include')
    if os.path.exists('/usr/local/opt/openblas/lib'):
        extra_link_args.append('/usr/local/opt/openblas/lib')
else:
    # OSX clang does not (yet) support openmp, so don't add it to compile
    # args
    extra_compile_args.append('-fopenmp')
    extra_link_args.append('-fopenmp')

extra_compile_args.append('-I'+os.path.join('src'))

ffi.set_source(
    '_prox_tv',
    """
    #include "TVopt.h"
    """,
    sources=sources,
    source_extension='.cpp',
    define_macros=[('NOMATLAB', 1)],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    libraries=['blas', 'lapack'],
    include_dirs=['/usr/include']
)


if __name__ == '__main__':
    ffi.compile(verbose=True)

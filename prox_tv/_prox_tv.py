r"""Proximal operators for 1d and 2d total variation problems.

The functions have the form ``tv<a>[w]_<b>d``, where:

    * ``<a>`` : 1 or 2, if the norm is :math:`\ell_1` or :math:`\ell_2`.
    * ``[w]`` : if present, the methods accepts a weighted norm.
    * ``<b>`` : 1 or 2, if the expected signals is 1- or 2-dimensional.

As the underlying library uses FORTRAN-style matrices (column-order), the given
matrices will be converted to this format if necessary.
"""

import numpy as np
import os
import os.path

from cffi import FFI

# The maximum number of returned info parameters.
_N_INFO = 3


_ffi = FFI()
_ffi.cdef("""
    typedef struct {
        ...;
    } Workspace;

    // Condat's implementatino.
    void TV1D_denoise(double* input, double* output, const int width,
                      const double lambda);
    // Ryan's implementation of Johnson's algorithm
    void dp(int n, double *y, double lam, double *beta);

    /* TV-L1 solvers */

    int tautString_TV1(double *y, double lambda, double *x,int n);
    int PN_TV1(double *y, double lambda, double *x, double *info, int n,
               double sigma, Workspace *ws);

    /* Weighted TV-L1 solvers */
    int PN_TV1_Weighted(double* Y, double* W, double* X, double* info, int n,
                        double sigma, Workspace* ws);
    int tautString_TV1_Weighted(double *y, double* lambda, double *x, int n);

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


# NOTE: The following hack is necessary, see
#       https://groups.google.com/forum/#!topic/python-cffi/puBLmTHBVmA
_cur_dir = os.getcwd()
os.chdir(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
_sources = ['src/' + fname for fname in (
    'condat_fast_tv.cpp', 'johnsonRyanTV.cpp', 'LPopt.cpp', 'TV2Dopt.cpp',
    'TV2DWopt.cpp', 'TVgenopt.cpp', 'TVL1opt.cpp', 'TVL1Wopt.cpp',
    'TVL2opt.cpp', 'TVLPopt.cpp', 'TVNDopt.cpp', 'utils.cpp'
)]

_lib = _ffi.verify(
    """
    typedef struct {
        /* Size of memory vectors */
        int n;
        /* Generic memory which can be used by 1D algorithms */
        double **d;
        int maxd, nd;
        int **i;
        int maxi, ni;
        /* Memory for inputs and outputs */
        double *in,*out;
        /* Warm restart variables */
        short warm;
        double *warmDual;
        double warmLambda;
    } Workspace;
    """,
    sources=_sources,
    define_macros=[('NOMATLAB', 1)],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
    libraries=['lapack'])

os.chdir(_cur_dir)


def _call(fn, *args):
    args_m = []
    for arg in args:
        if isinstance(arg, np.ndarray):
            args_m.append(_ffi.cast("double *", arg.ctypes.data))
        else:
            args_m.append(arg)
    return fn(*args_m)


def tv1_1d(x, w, sigma=0.05, method='tautstring'):
    r"""1D proximal operator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \sum_i |y_i - y_{i+1}|.

    Parameters
    ----------
    y : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    method : str
        The algorithm to be used, one of:

        * ``'tautstring'``
        * ``'pn'`` - projected Newton.
        * ``'condat'`` - Condat's segment construction method.
        * ``'dp'`` - Johnson's dynamic programming algorithm.

    sigma : float
        Tolerance for sufficient descent (used only if ``method='pn'``).

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert method in ('tautstring', 'pn', 'condat', 'dp')
    assert w >= 0
    y = np.zeros(np.size(x))
    if method == 'tautstring':
        _call(_lib.tautString_TV1, x, w, y, np.size(x))
    elif method == 'pn':
        info = np.zeros(_N_INFO)  # Holds [num of iterations, gap]
        _call(_lib.PN_TV1, x, w, y, info, np.size(x), sigma, _ffi.NULL)
    elif method == 'condat':
        _call(_lib.TV1D_denoise, x, y, np.size(x), w)
    else:
        _call(_lib.dp, np.size(x), x, w, y)
    return y


def tv1w_1d(x, w, method='tautstring', sigma=0.05):
    r"""Weighted 1D proximal operator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + \sum_i w_i |y_i - y_{i+1}|.

    Parameters
    ----------
    y : numpy array
        The signal we are approximating.
    w : numpy array
        The non-negative weights in the optimization problem.
    method : str
        Either ``'tautstring'`` or ``'pn'`` (projected Newton).
    sigma : float
        Tolerance for sufficient descent (used only if ``method='pn'``).

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert np.all(w >= 0)
    assert np.size(x) == np.size(w)
    y = np.zeros(np.size(x))
    if method == 'tautstring':
        _call(_lib.tautString_TV1_Weighted, x, w, y, np.size(x))
    else:
        info = np.zeros(_N_INFO)  # Holds [num of iterations, gap]
        _call(_lib.PN_TV1_Weighted,
              x, w, y, info, np.size(x), sigma, _ffi.NULL)
    return y


def tv2_1d(x, w, method='ms'):
    r"""1D proximal operator for :math:`\ell_2`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \sum_i (y_i - y_{i+1})^2.

    Parameters
    ----------
    y : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    method : str
        One of the following:

        * ``'ms'`` - More-Sorenson.
        * ``'pg'`` - Projected gradient.
        * ``'mspg'`` - More-Sorenson + projected gradient.

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert w >= 0
    assert method in ('ms', 'pg', 'mspg')
    info = np.zeros(_N_INFO)
    y = np.zeros(np.size(x), order='F')
    if method == 'ms':
        _call(_lib.more_TV2, x, w, y, info, np.size(x))
    elif method == 'pg':
        _call(_lib.PG_TV2, x, w, y, info, np.size(x))
    elif method == 'mspg':
        info = np.zeros(_N_INFO)
        _call(_lib.morePG_TV2, x, w, y, info, np.size(x), _ffi.NULL)
    return y


def tvp_1d(x, w, p, method='gp', max_iters=1000):
    r"""1D proximal operator for any :math:`\ell_p` norm.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \|y_i - y_{i+1}\|_p.

    Parameters
    ----------
    y : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    method : str
        The method to be used, one of the following:

         * ``'gp'`` - gradient projection
         * ``'ogp'`` - optimal gradient projection
         * ``'fista'`` - use the FISTA algorithm
         * ``'fw'`` - Frank-Wolfe
         * ``'gpfw'`` - hybrid gradient projection + Frank-Wolfe

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    methods = {
        'gp': _lib.GP_TVp,
        'ogp': _lib.OGP_TVp,
        'fista': _lib.FISTA_TVp,
        'fw': _lib.FW_TVp,
        'gpfw': _lib.GPFW_TVp,
    }
    assert method in methods
    assert w >= 0
    assert p >= 1
    info = np.zeros(_N_INFO)
    y = np.zeros(np.size(x), order='F')
    _call(methods[method], x, w, y, info, np.size(x), p, _ffi.NULL)
    return y


def tv1_2d(x, w, n_threads=8, max_iters=1000, method='yang'):
    r"""2D proximal oprator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 +
                       w \sum_{i,j} (|y_{i, j} - y_{i, j+1}| +
                       |y_{i,j} - y_{i+1,j}|).

    Parameters
    ----------
    y : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    str : method
        One of the following:

        * ``'yang'`` - Yang's algorithm.
        * ``'condat'`` - Condat's gradient.
        * ``'chambolle-pock'`` - Chambolle-Pock's gradient.

    n_threads : int
        Number of threads, used only for Douglas-Rachford.

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert w >= 0
    assert method in ('yang', 'condat', 'chambolle-pock')
    x = np.asfortranarray(x)
    y = np.asfortranarray(np.zeros(x.shape))
    info = np.zeros(_N_INFO)
    if method == 'yang':
        _call(_lib.Yang2_TV, x.shape[0], x.shape[1], x, w, y, max_iters, info)
    else:
        algorithm = 0 if method == 'condat' else 1
        _call(_lib.CondatChambollePock2_TV,
              x.shape[0], x.shape[1], x, w, y, algorithm, max_iters, info)

    return y


def tv1w_2d(x, w_col, w_row, max_iters=1000, n_threads=8):
    r"""2D weighted proximal operator for :math:`\ell_1` using DR splitting.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 +
                       \sum_{i,j} w^c_{i, j} |y_{i,j} - y_{i, j+1}| +
                       w^r_{i, j} |y_{i,j}-y_{i+1,j}|.

    Parameters
    ----------
    y : numpy array
        The MxN matrix we are approximating.
    w_col : numpy array
        The (M-1) x N matrix of column weights :math:`w^c`.
    w_row : numpy array
        The M x (N-1) matrix of row weights :math:`w^r`.

    Returns
    -------
    numpy array
        The MxN solution of the optimization problem.
    """
    assert np.all(w_col >= 0)
    assert np.all(w_row >= 0)
    M, N = x.shape
    assert w_col.shape == (M-1, N)
    assert w_row.shape == (M, N-1)
    x = np.asfortranarray(x)
    y = np.zeros(x.shape, order='F')
    w_col = np.asfortranarray(w_col)
    w_row = np.asfortranarray(w_row)
    info = np.zeros(_N_INFO)
    _call(_lib.DR2L1W_TV, M, N, x, w_col, w_row, y, n_threads, max_iters, info)
    return y


def tvp_2d(x, w_col, w_row, p_col, p_row, n_threads=8, max_iters=1000):
    r"""2D proximal operator for any :math:`\ell_p` norm.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2}\|x-y\|^2 + w^r \|D_\mathrm{row}(y)\|_{p_1} +
                                              w^c \|D_\mathrm{col}(y) \|_{p_2},

    where :math:`\mathrm D_{row}` and :math:`\mathrm D_{col}` take the
    differences accross rows and columns respectively.

    Parameters
    ----------
    y : numpy array
        The matrix signal we are approximating.
    p_col : float
        Column norm.
    p_row : float
        Row norm.
    w_col : float
        Column penalty.
    w_row : float
        Row penalty.

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert w_col >= 0
    assert w_row >= 0
    assert p_col >= 1
    assert p_row >= 1

    info = np.zeros(_N_INFO)
    x = np.asfortranarray(x)
    y = np.zeros(np.shape(x), order='F')
    _call(_lib.DR2_TV,
          x.shape[0], x.shape[1], x, w_col, w_row, p_col, p_row, y,
          n_threads, max_iters, info)
    return y

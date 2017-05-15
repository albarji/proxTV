r"""**proxTV** is a toolbox implementing blazing fast implementations of Total Variation proximity operators.

The library provides efficient solvers for the following Total Variation proximity problems:

**Standard (l1) Total Variation on a 1-dimensional signal**

    .. figure::  img/TV1.png

**Quadratic (l2) Total Variation on a 1-dimensional signal**

    .. figure::  img/TV2.png

**lp-norm Total Variation on a 1-dimensional signal**

    .. figure::  img/TVp.png

**Weighted Total Variation on a 1-dimensional signal**

    .. figure::  img/TV1w.png

**Anisotropic Total Variation on a 2-dimensional signal**

    .. figure::  img/TV2D.png

**lp-norm Anisotropic Total Variation on a 2-dimensional signal**

    .. figure::  img/TV2Dp.png

**Weighted Anisotropic Total Variation on a 2-dimensional signal**

    .. figure::  img/TV2Dw.png

**Anisotropic Total Variation on a 3-dimensional signal**

    .. figure::  img/TV3D.png

**Generalized N-dimensional Anisotropic Total Variation**

    .. figure::  img/TVND.png
    (with X(di) every possible 1-dimensional slice of X following dimension di)


As a general rule the functions in this package have the form ``tv<a>[w]_<b>d``, where:

    * ``<a>`` : 1, 2 or p, if the norm is :math:`\ell_1`, :math:`\ell_2` or the general :math:`\ell_p`
    * ``[w]`` : if present, the methods accepts a weighted norm.
    * ``<b>`` : 1 or 2, if the expected signals is 1- or 2-dimensional.

The only expection is the function **tvgen** that solves generalized
Total Variation problems, recommended only to advanced users.

As the underlying library uses FORTRAN-style matrices (column-order), the given
matrices will be converted to this format if necessary.

If you find this toolbox useful please reference the following papers:

* Fast Newton-type Methods for Total Variation Regularization. Alvaro Barbero, Suvrit Sra. ICML 2011 proceedings.

* Modular proximal optimization for multidimensional total-variation regularization. Alvaro Barbero, Suvrit Sra. http://arxiv.org/abs/1411.0589

"""

import numpy as np
from _prox_tv import ffi, lib

# The maximum number of returned info parameters.
_N_INFO = 3


def _call(fn, *args):
    args_m = []
    for arg in args:
        if isinstance(arg, np.ndarray):
            args_m.append(ffi.cast("double *", arg.ctypes.data))
        else:
            args_m.append(arg)
    return fn(*args_m)


def force_float_scalar(x):
    r"""Forces an scalar value into float format

    Parameters
    ----------
    x: scalar value to check

    Returns
    -------
    float
        Float representation of the provided value
    """
    if not isinstance(x, float):
        return float(x)
    else:
        return x


def force_float_matrix(x):
    r"""Forces a numpy matrix into float format

    Parameters
    ----------
    x: numpy array
        matrix to check

    Returns
    -------
    numpy array
        Float representation of the provided matrix
    """
    # Check input is numpy array
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except Exception:
            raise TypeError("Input must be a numpy matrix or compatible object")
    # Enforce float type
    if x.dtype != np.dtype('float64'):
        return x.astype('float')
    else:
        return x


def tv1_1d(x, w, method='hybridtautstring', sigma=0.05, maxbacktracks=None):
    r"""1D proximal operator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \sum_i |y_i - y_{i+1}|.

    Parameters
    ----------
    x : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    method : str
        The algorithm to be used, one of:

        * ``'classictautstring'`` - classic Taut String method
        * ``'linearizedtautstring'`` - linearized Taut String method
        * ``'hybridtautstring'`` - hybrid classic+linear Taut String method
        * ``'pn'`` - projected Newton.
        * ``'condat'`` - Condat's segment construction method.
        * ``'dp'`` - Johnson's dynamic programming algorithm.
        * ``'condattautstring'`` - Condat's implementation of classic taut
            string method.
        * ``'kolgomorov'`` - Kolmogorov et al's message passing method.

    sigma : float
        Tolerance for sufficient descent (used only if ``method='pn'``).
        
    maxbacktracks: float
        Backtrack steps before switching (used only if ``method='hybridtautstring'``)

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    METHODS = {
        'classictautstring' : _classictautstring_TV1,
        'linearizedtautstring' : _linearizedtautstring_TV1,
        'hybridtautstring' : _hybridtautstring_TV1,
        'pn' : _pn_TV1,
        'condat' : _condat_TV1,
        'dp' : _dp_TV1,
        'condattautstring' : _condattautstring_TV1,
        'kolmogorov' : _kolmogorov_TV1
    }
    assert method in METHODS
    assert w >= 0
    w = force_float_scalar(w)
    x = force_float_matrix(x)
    y = np.zeros(np.size(x))
    METHODS[method](x, w, y, sigma=sigma, maxbacktracks=maxbacktracks)
    return y

def _classictautstring_TV1(x, w, y, **kwargs):
    """Classic taut string method for TV1 proximity"""
    _call(lib.classicTautString_TV1, x, np.size(x), w, y)

def _linearizedtautstring_TV1(x, w, y, **kwargs):
    """Linearized taut string method for TV1 proximity"""
    _call(lib.linearizedTautString_TV1, x, w, y, np.size(x))

def _hybridtautstring_TV1(x, w, y, **kwargs):
    """Hybrid taut string method for TV1 proximity"""
    if kwargs['maxbacktracks'] is None:
        _call(lib.hybridTautString_TV1, x, np.size(x), w, y)    
    else:
        _call(lib.hybridTautString_TV1_custom, x, np.size(x), w, y, 
              kwargs['maxbacktracks']) 

def _pn_TV1(x, w, y, **kwargs):
    """Projected Newton method for TV1 proximity"""
    info = np.zeros(_N_INFO)  # Holds [num of iterations, gap]
    _call(lib.PN_TV1, x, w, y, info, np.size(x), kwargs['sigma'], ffi.NULL)

def _condat_TV1(x, w, y, **kwargs):
    """Condat's method for TV1 proximity"""
    _call(lib.TV1D_denoise, x, y, np.size(x), w)

def _condattautstring_TV1(x, w, y, **kwargs):
    """Condat's implementation of the taut string method for TV1 proximity"""
    _call(lib.TV1D_denoise_tautstring, x, y, np.size(x), w)

def _kolmogorov_TV1(x, w, y, **kwargs):
    """Kolmogorov's method for TV1 proximity"""
    _call(lib.SolveTVConvexQuadratic_a1_nw, np.size(x), x, w, y)

def _dp_TV1(x, w, y, **kwargs):
    """Johnson's dynamic programming method for TV1 proximity"""
    _call(lib.dp, np.size(x), x, w, y)

def tv1w_1d(x, w, method='tautstring', sigma=0.05):
    r"""Weighted 1D proximal operator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + \sum_i w_i |y_i - y_{i+1}|.

    Parameters
    ----------
    x : numpy array
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
    assert np.size(x)-1 == np.size(w)
    w = force_float_matrix(w)
    x = force_float_matrix(x)
    y = np.zeros(np.size(x))
    if method == 'tautstring':
        _call(lib.tautString_TV1_Weighted, x, w, y, np.size(x))
    else:
        info = np.zeros(_N_INFO)  # Holds [num of iterations, gap]
        _call(lib.PN_TV1_Weighted,
              x, w, y, info, np.size(x), sigma, ffi.NULL)
    return y


def tv2_1d(x, w, method='mspg'):
    r"""1D proximal operator for :math:`\ell_2`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \sum_i (y_i - y_{i+1})^2.

    Parameters
    ----------
    x : numpy array
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
    METHODS = {
        'ms' : _ms_tv2,
        'pg' : _pg_tv2,
        'mspg': _mspg_tv2
    }
    assert w >= 0
    assert method in METHODS

    w = force_float_scalar(w)
    x = force_float_matrix(x)
    info = np.zeros(_N_INFO)
    y = np.zeros(np.size(x), order='F')
    METHODS[method](x, w, y, info)
    return y

def _ms_tv2(x, w, y, info):
    """More-Sorensen method for TV2 proximity"""
    _call(lib.more_TV2, x, w, y, info, np.size(x))

def _pg_tv2(x, w, y, info):
    """Projected Gradient method for TV2 proximity"""
    _call(lib.PG_TV2, x, w, y, info, np.size(x))

def _mspg_tv2(x, w, y, info):
    """More-Sorensen + Projected Gradient hybrid method for TV2 proximity"""
    _call(lib.morePG_TV2, x, w, y, info, np.size(x), ffi.NULL)

def tvp_1d(x, w, p, method='gpfw', max_iters=0):
    r"""1D proximal operator for any :math:`\ell_p` norm.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 + w \|y_i - y_{i+1}\|_p.

    Parameters
    ----------
    x : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    method : str
        The method to be used, one of the following:

         * ``'gp'`` - gradient projection
         * ``'fw'`` - Frank-Wolfe
         * ``'gpfw'`` - hybrid gradient projection + Frank-Wolfe

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    methods = {
        'gp': lib.GP_TVp,
        'fw': lib.FW_TVp,
        'gpfw': lib.GPFW_TVp,
    }
    assert method in methods
    assert w >= 0
    assert p >= 1
    w = force_float_scalar(w)
    p = force_float_scalar(p)
    x = force_float_matrix(x)
    info = np.zeros(_N_INFO)
    y = np.zeros(np.size(x), order='F')
    _call(methods[method], x, w, y, info, np.size(x), p, ffi.NULL)
    return y


def tv1_2d(x, w, n_threads=1, max_iters=0, method='dr'):
    r"""2D proximal oprator for :math:`\ell_1`.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 +
                       w \sum_{i,j} (|y_{i, j} - y_{i, j+1}| +
                       |y_{i,j} - y_{i+1,j}|).

    Parameters
    ----------
    x : numpy array
        The signal we are approximating.
    w : float
        The non-negative weight in the optimization problem.
    str : method
        One of the following:

        * ``'dr'`` - Douglas Rachford splitting.
        * ``'pd'`` - Proximal Dykstra splitting.
        * ``'yang'`` - Yang's algorithm.
        * ``'condat'`` - Condat's gradient.
        * ``'chambolle-pock'`` - Chambolle-Pock's gradient.
        * ``'kolmogorov'`` - Kolmogorov's splitting.

    n_threads : int
        Number of threads, used only for Proximal Dykstra
            and Douglas-Rachford.

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    METHODS = {
        'yang' : _yang_tv2d,
        'dr' : _dr_tv2d,
        'pd' : _pd_tv2d,
        'kolmogorov' : _kolmogorov_tv2d,
        'condat': _condat_tv2d,
        'chambolle-pock':  _chambollepock_tv2d,
        'chambolle-pock-acc': _chambollepockacc_tv2d,
    }
    assert w >= 0
    assert method in METHODS
    x = np.asfortranarray(x, dtype='float64')
    w = force_float_scalar(w)
    y = np.asfortranarray(np.zeros(x.shape))
    info = np.zeros(_N_INFO)
    METHODS[method](x, w, y, max_iters, info, n_threads=n_threads)
    return y

def _yang_tv2d(x, w, y, max_iters, info, **kwargs):
    """Yang's method for 2D TV proximity"""
    _call(lib.Yang2_TV, x.shape[0], x.shape[1], x, w, y, max_iters, info)

def _dr_tv2d(x, w, y, max_iters, info, **kwargs):
    """Douglas Rachford parallel projections method for 2D TV proximity"""
    _call(lib.DR2_TV, x.shape[0], x.shape[1], x, w, w, 1.0, 1.0, y,
              kwargs['n_threads'], max_iters, info)

def _pd_tv2d(x, w, y, max_iters, info, **kwargs):
    """Proximal dykstra method for 2D TV proximity"""
    _call(lib.PD2_TV, x, (w, w), (1, 1), (1, 2), y, info, x.shape, 2, 2,
              kwargs['n_threads'], max_iters)

def _kolmogorov_tv2d(x, w, y, max_iters, info, **kwargs):
    """Kolmogorov's method for 2D TV proximity"""
    _call(lib.Kolmogorov2_TV, x.shape[0], x.shape[1], x, w, y, max_iters, 
              info)

def _condat_tv2d(x, w, y, max_iters, info, **kwargs):
    """Condat's method for 2D TV proximity"""
    _condatchambollepock_tv2d(x, w, y, max_iters, info, algorithm=0)

def _chambollepock_tv2d(x, w, y, max_iters, info, **kwargs):
    """Chambolle and Pock's method for 2D TV proximity"""
    _condatchambollepock_tv2d(x, w, y, max_iters, info, algorithm=1)

def _chambollepockacc_tv2d(x, w, y, max_iters, info, **kwargs):
    """Accelerated Chambolle and Pock's method for 2D TV proximity"""
    _condatchambollepock_tv2d(x, w, y, max_iters, info, algorithm=2)

def _condatchambollepock_tv2d(x, w, y, max_iters, info, **kwargs):
    """Wrapper function for general implementaton of Chambolle-Pock + Condat"""
    _call(lib.CondatChambollePock2_TV, x.shape[0], x.shape[1], x, w, y, 
          kwargs['algorithm'], max_iters, info)

def tv1w_2d(x, w_col, w_row, max_iters=0, n_threads=1):
    r"""2D weighted proximal operator for :math:`\ell_1` using DR splitting.

    Specifically, this optimizes the following program:

    .. math::

        \mathrm{min}_y \frac{1}{2} \|x-y\|^2 +
                       \sum_{i,j} w^c_{i, j} |y_{i,j} - y_{i, j+1}| +
                       w^r_{i, j} |y_{i,j}-y_{i+1,j}|.

    Parameters
    ----------
    x : numpy array
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
    x = np.asfortranarray(x, dtype='float64')
    y = np.zeros(x.shape, order='F')
    w_col = np.asfortranarray(w_col, dtype='float64')
    w_row = np.asfortranarray(w_row, dtype='float64')
    info = np.zeros(_N_INFO)
    _call(lib.DR2L1W_TV, M, N, x, w_col, w_row, y, n_threads, max_iters, info)
    return y


def tvp_2d(x, w_col, w_row, p_col, p_row, n_threads=1, max_iters=0):
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
    x = np.asfortranarray(x, dtype='float64')
    w_col = force_float_scalar(w_col)
    w_row = force_float_scalar(w_row)
    p_col = force_float_scalar(p_col)
    p_row = force_float_scalar(p_row)
    y = np.zeros(np.shape(x), order='F')
    _call(lib.DR2_TV,
          x.shape[0], x.shape[1], x, w_col, w_row, p_col, p_row, y,
          n_threads, max_iters, info)
    return y


def tvgen(x, ws, ds, ps, n_threads=1, max_iters=0):
    r"""General TV proximal operator for multidimensional signals

    Specifically, this optimizes the following program:

    .. math::

      \min_X \frac{1}{2} ||X-Y||^2_2 + \sum_i w_i \sum_j TV^{1D}(X(d_i)_j,p_i)

    where :math:`X(d_i)_j` every possible 1-dimensional fiber of X following
    the dimension d_i, and :math:`TV^{1D}(z,p)` the 1-dimensional
    :math:`\ell_p`-norm total variation over the given fiber :math:`z`.
    The user can specify the number :math:`i` of penalty terms.

    Parameters
    ----------
    x : numpy array
        The matrix signal we are approximating.
    ws : iterable
        Weights to apply in each penalty term.
    ds : iterable
        Dimensions over which to apply each penalty term.
        Must be of equal length to ws.
    ps : iterable
        Norms to apply in each penalty term.
        Must be of equal length to ws.
    n_threads : int
        number of threads to use in the computation
    max_iters : int
        maximum number of iterations to run the solver.

    Returns
    -------
    numpy array
        The solution of the optimization problem.
    """
    assert len(ws) == len(ds)
    assert len(ws) == len(ps)
    assert n_threads >= 1
    assert max_iters >= 0

    info = np.zeros(_N_INFO)
    x = np.asfortranarray(x, dtype='float64')
    ws = force_float_matrix(ws)
    ps = force_float_matrix(ps)
    y = np.zeros(np.shape(x), order='F')

    # Run algorithm depending on the structure of the data and the requested
    # penalties.

    # Bidimensional signal with one penalty term across each dimension (full
    # 2-dimensional TV proximity): Douglas-Rachford splitting
    if len(ds) == 2 & ds[0] == 1 & ds[1] == 2:
        _call(lib.DR2_TV,
              x.shape[0], x.shape[1], x, ws[0], ws[1], ps[0], ps[1], y,
              n_threads, max_iters, info)
    # 2 arbitrary terms: Proximal Dykstra
    elif len(ws) == 2:
        _call(lib.PD2_TV,
              x, ws, ps, ds, y, info, x.shape, len(x.shape), 2, n_threads,
              max_iters)
    # More terms: Parallel Proximal Dykstra
    else:
        _call(lib.PD_TV,
              x, ws, ps, ds, y, info, x.shape, len(x.shape), len(ws),
              n_threads, max_iters)

    return y

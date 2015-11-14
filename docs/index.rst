Proximal total-variation operators
**********************************

.. automodule:: prox_tv

Function reference
==================

One dimensional total variation problems
----------------------------------------


.. autofunction:: prox_tv.tv1_1d

.. autofunction:: prox_tv.tv1w_1d

.. autofunction:: prox_tv.tv2_1d

.. autofunction:: prox_tv.tvp_1d


Two dimensional total variation problems
----------------------------------------


.. autofunction:: prox_tv.tv1_2d

.. autofunction:: prox_tv.tv1w_2d

.. autofunction:: prox_tv.tvp_2d

Generalized total variation problems
------------------------------------

.. autofunction:: prox_tv.tvgen

Examples
========

1D examples
-----------

Filter 1D signal using TV-L1 norm::

    tv1_1d(x, w)
    
Filter 1D signal using weighted TV-L1 norm (for x vector of length N, weights vector of length N-1)::

    tv1w_1d(x, weights)
    
Filter 1D signal using TV-L2 norm::

    tv2_1d(x, w)
    
Filter 1D signal using both TV-L1 and TV-L2 norms::

    tvgen(X, [w1, w2], [1, 1], [1, 2])
    
2D examples
-----------

Filter 2D signal using TV-L1 norm::

    tv1_2d(X, w)
    
or::
        
    tvgen(X, [w, w], [1, 2], [1, 1])

Filter 2D signal using TV-L2 norm::

    tvp_2d(X, w)
    
or::

    tvgen(X, [w, w], [1, 2], [2, 2])
    
Filter 2D signal using 4 parallel threads::

    tv1_2d(X, w, n_threads=4)
    
or::

    tvgen(X, [w, w], [1, 2], [1, 1], n_threads=4)

Filter 2D signal using TV-L1 norm for the rows, TV-L2 for the columns, and different penalties::

    tvgen(X, [wRows, wCols], [1, 2], [1, 2])
    
Filter 2D signal using both TV-L1 and TV-L2 norms::

    tvgen(X, [w1, w1, w2, w2], [1, 2, 1, 2], [1, 1, 2, 2])
    
Filter 2D signal using weighted TV-L1 norm (for X image of size MxN, W1 weights of size (M-1)xN, W2 weights of size Mx(N-1))::

    tv1w_2d(X, W1, W2)
    
3D examples
-----------

Filter 3D signal using TV-L1 norm::

    tvgen(X, [w, w, w], [1, 2, 3], [1, 1, 1])

Filter 3D signal using TV-L2 norm, not penalizing over the second dimension::

    tvgen(X , [w, w], [1, 3], [2, 2])
    

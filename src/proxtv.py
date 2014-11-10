#
# Python module for using proxTV library functions
#
# (c) Copyright 2014  Alvaro Barbero and Suvrit Sra
#

import proxtv_internal as tv
import numpy as np

#
# [x,info] = TV(y,lambda,p,threads)  Proximity operator for Lp Total Variation.
#
# Solves the proximity problem associated with the n-dimensional Total Variation Lp norm.
# Depending on the dimension and norm of choice, a different algorithm is used for the
# optimization.
# Any norm p>=1 is accepted.
#
# Inputs:
#   - y: numpy multi-array, input of the proximity operator.
#   - lambda: premultiplier of the norm.
#   - [p]: norm (default 1).
#   - [threads]: number of threads (default 1). Used only for 2-D or larger signals.
#   - [mit]: max num iters for proximal Dykstra (default 10).
#
# Outputs:
#   - x: numpy multi-array, solution of the proximity problem.
#   - info: statistical info of the run algorithm (if any)
#       info.iters: number of iterations run (major iterations for the 2D case)
#       info.stop: value of the stopping criterion.
def TV(X, lambda, p, nthreads, mit)
    # Check input is numpy array
    if type(X) != np.ndarray:
        raise ProxTVInputArgumentException('X','Input signal must be numpy.ndarray object')
    
    # Check dimension of input signal
    dim = X.ndim
    
    #TODO
        
#
# Classes for proxTV exceptions
#

class ProxTVException(Exception):
    """Base class for exceptions in proxTV"""
    pass
    
class ProxTVInputArgumentException(ProxTVException):
    """Exception raised for errors the inputs when invoking proxTV functions.
    
    Attributes:
        expr -- input expression in which the error occurred
        msg -- explanation of the error
    """
    
    def __init__(self, expr, msg)
        self.expr = expr
        self.msg = msg


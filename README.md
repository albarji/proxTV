proxTV
======

Matlab toolbox for fast Total Variation proximity operators.

For an up-to-date version, check https://github.com/albarji/proxTV .

Index
-----

1. Quick start guide.
2. Referencing.
3. Installation.
4. Usage.
5. Examples.
6. Demos.
7. Contact.
8. Acknowledgements

1. Quick start guide
--------------------

To install proxTV just type "install" at the Matlab prompt once located at proxTV folder. If any problem arises please refer to the "Installation" section in this file.

After that the TV solver can be invoked easily through the general purpose "TV" function. For instance,

    TV(X,lambda)
    
solves TV for X signal and lambda regularizer. The dimensionality of X is automatically checked and an adequate solver is applied.

To solve TV for a general TV-Lp norm just add the value of p as a third argument,

    TV(X,lambda,p)
    
Weighted versions of TV can also be solved by using exactly the same interface, but providing a vector of lambda weights instead of a scalar. For multidimensional signals the relevant weights are provided as a cell array; the "Usage" section for more detailts on this and more advanced uses of toolbox.

2. Referencing
--------------

If you find this toolbox useful please reference the following papers:

    Fast Newton-type Methods for Total Variation Regularization. Álvaro Barbero, Suvrit Sra. ICML 2011 proceedings.
    Modular proximal optimization for multidimensional total-variation regularization. Álvaro Barbero, Suvrit Sra. Online.
    
whose Bibtex entries are

    @inproceedings{conf/icml/Barbero11,
      added-at = {2011-12-16T00:00:00.000+0100},
      author = {Barbero, Álvaro and Sra, Suvrit},
      biburl = {http://www.bibsonomy.org/bibtex/214ce9f5c15d1d462bd264d8af9e4c3c7/dblp},
      booktitle = {ICML},
      crossref = {conf/icml/2011},
      editor = {Getoor, Lise and Scheffer, Tobias},
      interhash = {5d6359b6c7f4d0fb6de36aada6827a3e},
      intrahash = {14ce9f5c15d1d462bd264d8af9e4c3c7},
      keywords = {dblp},
      pages = {313-320},
      publisher = {Omnipress},
      timestamp = {2011-12-16T00:00:00.000+0100},
      title = {Fast Newton-type Methods for Total Variation Regularization.},
      url = {http://dblp.uni-trier.de/db/conf/icml/icml2011.html#JimenezS11},
      year = 2011
    }
    @Article{barberoTV14,
      Title                    = {Modular proximal optimization for multidimensional total-variation regularization},
      Author                   = {\'Alvaro Barbero and Suvrit Sra},
      Year                     = {2014},
      Url                      = {http://arxiv.org/abs/1411.0589}
    }

3. Installation
---------------

To install proxTV follow the steps:

1) Open Matlab.
2) Change directory to the folder where this README file is located.
3) Type: "install"
4) The message "proxTV successfully installed" will appear. If instead an error message shows up, check your mex compiler configuration (type "mex -setup"). If the error was produced by the "openmp" library, you might need to install it, or you can install proxTV with no multi-thread features by typing "install 0". 
5) The install script automatically adds the relevant proxTV folders to your Matlab path. If this were to fail for some reason, you should manually add the /src subfolder to your Matlab path. You can do this e.g. by using the "pathtool()" utility. You may also add the /demos subfolder to try out the included demos.

Note: this toolbox has only been tested under Linux. Installation might require LAPACK (http://www.netlib.org/lapack/) and BLAS (http://www.netlib.org/blas/) libraries.

4. Usage
--------

Two main functions conform the proxTV toolbox: TV and TVgen. The first one provides basic options over the Total Variation problem, while the second one allows a more advanced configuration. In general, the TV function should suffice for most uses.

a) TV
·····

Solves Total Variation proximity operators for n-dimensional signals, applying a TV-Lp norm. The inputs and outputs of this function are:

    [x,info] = TV(y,lambda,p,threads)

 Inputs:
   - y: input of the proximity operator.
   - lambda: premultiplier of the norm.
   - (Optional) p: norm. Default is p = 1.
   - (Optional) threads: number of threads (default 1). Used only for 2-D or higher-dimensional signals.

 Outputs:
   - x: solution of the proximity problem.
   - info: statistical info of the run algorithm:
       info.iters: number of iterations run (major iterations for the 2D case)
       info.stop: value of the stopping criterion.
       
For 1-dimensional signals the problem solved is

    min_x  1/2 ||x-y||_2^2 + lambda * TV(x,p)
    
  = min_x  1/2 ||x-y||_2^2 + lambda * ||sum_i (x_i - x_(i-1))||_p
    
Using p=1 results in the standard Total Variation regularizer, which generates a "blocky" reconstruction of the signal. Using p=2 instead produces a smoother reconstruction.

For 2-dimensional signals (e.g. images) the problem solved is

    min_X  1/2 ||X(:)-Y(:)||_2^2 + lambda * ( sum_i TV(Xri,p) + sum_j TV(Xcj,p) )
    
where X(:) is the vectorization (either along the rows or the columns) of X, Xri is the i-th row of X and Xcj is the j-th column of X. Using p=1 results in an anisotropic denoising filter.

For D-dimensional signals the problem being solved becomes

    min_X 0.5 ||X(:)-Y(:)||^2 + lambda * ( sum_i1 TV(X[1,i1],p) + sum_i2 TV(X[2,i2],p) + ... + sum_iD TV(X[D,iD],D) )
    
where X[i,d] is the i-th 1-dimensional fiber of X along its d-th dimension.

If a vector of weights is provided for the lambda parameter instead of an scalar value, the special weighted version of TV is solved,

    min_x 0.5 ||x-y||^2 + sum_i lambda_i |x_i - x_(i-1)|
    
were each difference among signal entries x_i and x_(i-1) is penalized using a different weight lambda_i. In the case of multidimensional signals the weighted problem is

    min_X 0.5 ||X(:)-Y(:)||^2 + sum_i1 TV(X[1,i1],p,W1[i1]) + sum_i2 TV(X[2,i2],p,W2_i) + ... + sum_iD TV(X[D,iD],D,WD_i) )
    
where Wd[i] is the 1-dimensional fiber of weights along the d-th dimension applied to X[i,d]. Weight tensors are provided in TV function as the lambda parameter through a cell array in the form {W1, W2, ..., Wd} (see the examples in the "Examples" section)

b) TVgen
········

Solves a generalized TV proximity operator for a multidimensional signal, in the form

    min_x 0.5 ||x-y||^2 + sum_i P(x,lambda_i,d_i,p_i)
    
where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional fiber of x following the dimension d_i. The inputs and outputs of this function are:

 Inputs:
   - y: input signal.
   - lambdas: vector of lambda penalties of each penalty term.
   - ds: vector of dimensions of application of each penalty term.
   - norms: vector of norms of each penalty term.
   - (Optional) threads: number of threads to use (default: 1)

 Outputs:
   - x: solution of the proximity problem.
   - info: statistical info of the run algorithm:
       info.iters: number of major iterations run.
       info.stop: value of the stopping criterion.
       
When possible, TV should be preferred. See the Examples section next for some specific examples on using this function.

5. Examples
-----------

1D examples
···········

- Filter 1D signal using TV-L1 norm:
    TV(x,lambda)
    
- Filter 1D signal using weighted TV-L1 norm (for x vector of length N, weights vector of length N-1)
    TV(x,weights)
    
- Filter 1D signal using TV-L2 norm:
    TV(x,lambda,2)
    
- Filter 1D signal using both TV-L1 and TV-L2 norms:
    TVgen(X,[lambda1 lambda2],[1 1],[1 2])
    
2D examples
··········· 

- Filter 2D signal using TV-L1 norm:
    TV(X,lambda)
        or
    TVgen(X,[lambda lambda],[1 2],[1 1])

- Filter 2D signal using TV-L2 norm:
    TV(X,lambda,2)
        or
    TVgen(X,[lambda lambda],[1 2],[2 2])
    
- Filter 2D signal using 4 parallel threads (last argument):
    TV(X,lambda,1,4)
        or
    TVgen(X,[lambda lambda],[1 2],[1 1],4)

- Filter 2D signal using TV-L1 norm for the rows, TV-L2 for the columns, and different penalties:
    TVgen(X,[lambdaRows lambdaCols],[1 2],[1 2])
    
- Filter 2D signal using both TV-L1 and TV-L2 norms:
    TVgen(X,[lambda1 lambda1 lambda2 lambda2],[1 2 1 2],[1 1 2 2])
    
- Filter 2D signal using weighted TV-L1 norm (for X image of size MxN, W1 weights of size (M-1)xN, W2 weights of size Mx(N-1))
    TV(X, {W1, W2})
    
3D examples
··········· 

- Filter 3D signal using TV-L1 norm:
    TV(X,lambda)
        or
    TVgen(X,[lambda lambda lambda],[1 2 3],[1 1 1])

- Filter 3D signal using TV-L2 norm, not penalizing over the second dimension:
    TVgen(X,[lambda lambda],[1 3],[2 2])

    
6. Demos
--------

Some demos in the form of Matlab scripts showing how to work with proxTV are included in the subfolder /demos. They are:

- demo_filter_signal: TV-L1, TV-L2 and weighted TV-L1 filtering of 1-dimensional signals.
- demo_filter_image: TV-L1 filtering of 2-dimensional image.
- demo_filter_image_color: TV-L1 filtering of 3-dimensional image (length, width and color).
- demo_filter_image_threads: multi-thread TV-L1 filtering of 2-dimensional image.

7. Contact
----------

For any questions and comments, please email alvaro.barbero@uam.es

8. Acknowledgements
-------------------

We wish to thank Zico Kolter for pointing out a bug in version 1.0 of this code.

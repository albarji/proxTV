% Given a reference multidimensional signal y and a series of penalty terms P(x,lambda,d,p), solves the generalized Total Variation
% proximity operator
%    
%        min_x 0.5 ||x-y||^2 + sum_i P(x,lambda_i,d_i,p_i) .
%        
% where P(x,lambda_i,d_i,p_i) = lambda_i * sum_j TV(x(d_i)_j,p_i) with x(d)_j every possible 1-dimensional slice of x
% following the dimension d_i, TV(x,p) the TV-Lp prox operator for x.
%
% Inputs:
%   - y: input signal.
%   - lambdas: vector of lambda penalties of each penalty term.
%   - ds: vector of dimensions of application of each penalty term.
%   - norms: vector of norms of each penalty term.
%   - [threads]: number of threads to use (default: 1)
%
% Outputs:
%   - x: solution of the proximity problem.
%   - info: statistical info of the run algorithm:
%       info.iters: number of major iterations run.
%       info.stop: value of the stopping criterion.
%
% Examples:
%
%   - Filter 2D signal using TV-L1 norm:
%        TVgen(X,[lambda lambda],[1 2],[1 1])
%
%   - Filter 2D signal using TV-L2 norm:
%        TVgen(X,[lambda lambda],[1 2],[2 2])
%
%   - Filter 2D signal using TV-L1 norm for the rows, TV-L2 for the columns, and different penalties:
%        TVgen(X,[lambdaRows lambdaCols],[1 2],[1 2])
%
%   - Filter 1D signal using both TV-L1 and TV-L2 norms:
%        TVgen(X,[lambda1 lambda2],[1 1],[1 2])
%
%   - Filter 3D signal using TV-L1 norm:
%        TVgen(X,[lambda lambda lambda],[1 2 3],[1 1 1])
%
%   - Filter 3D signal using TV-L2 norm, not penalizing over the second dimension:
%        TVgen(X,[lambda lambda],[1 3],[2 2])
%
%   - Filter 2D signal using both TV-L1 and TV-L2 norms:
%        TVgen(X,[lambda1 lambda1 lambda2 lambda2],[1 2 1 2],[1 1 2 2])
%
%   - ... and so on, any combination of norms and dimensions is possible.
function [x,info] = TVgen(y,lambdas,ds,norms,threads)
    % Check inputs
    if nargin < 5, threads = 1; end;
    if length(lambdas) ~= length(ds) || length(ds) ~= length(norms)
        fprintf(1,'ERROR (prox_TVgen): arguments defining penalties differ in length.\n');
        x = y;
        return;
    end;
    if sum(norms >= 1) ~= length(norms)
        fprintf(1,'ERROR (prox_TVgen): norms should be >=1 or inf.\n');
        x = y;
        return;
    end;

    % Invoke C solver
    [x,in] = solveTVgen(y,lambdas,ds,norms,threads);
    info.iters = in(1);
    info.stop = in(2);
end

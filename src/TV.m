%
% [x,info] = TV(y,lambda,p,threads)  Proximity operator for Lp Total Variation.
%
% Solves the proximity problem associated with the n-dimensional Total Variation Lp norm.
% Depending on the dimension and norm of choice, a different algorithm is used for the
% optimization.
% Any norm p>=1 is accepted.
%
% Inputs:
%   - y: input of the proximity operator.
%   - lambda: premultiplier of the norm.
%   - [p]: norm (default 1).
%   - [threads]: number of threads (default 1). Used only for 2-D or larger signals.
%   - [mit]: max num iters for combiner algorithm (default: as defined by the algorithm)).
%
% Outputs:
%   - x: solution of the proximity problem.
%   - info: statistical info of the run algorithm (if any)
%       info.iters: number of iterations run (major iterations for the 2D case)
%       info.stop: value of the stopping criterion.
%
function [x,info] = TV(y,lambda,p,threads,mit)
    % Check inputs
    if (~exist('p', 'var')), p = 1;  end;
    if (~exist('threads', 'var')), threads=1; end;
    if (~exist('mit', 'var')), mit=0; end; % 0 means let the algorithm use it's default value
    
    % Case where weights are given in a cell (weighted N-dimensional filtering case)
    if iscell(lambda)
        % Check compatible dimensionality with inputs
        if ndims(y) ~= length(lambda)
            error 'for an N-dimensional signal the weights must be provided in a cell array of length N'
        end
        % Currently only 2D weighted filtering is supported
        if length(lambda) ~= 2 
            error 'only 1D and 2D weighted filtering is supported'
        end
        % Currently only L1 norm is amenable for weighted TV
        if ( p ~= 1 )
            error 'Currently only L1 norm is accepted for weighted TV'
        end
        % Extract weights
        W1 = lambda{1};
        W2 = lambda{2};
        % Apply solver
        [x, in] = solveTV2DL1W(y, W1, W2, threads, mit);
        return;
    end
    
    % Special case where a vector of weights is provided: use weighted TV version
    if ( length(lambda) > 1 )
        % Check there are enough weights
        if ( length(y)-1 ~= length(lambda) )
            error 'lambda should be either a scalar or a weights vector with length(lambda) == length(y)-1'
        end
        % Currently only 1D signals are amenable for weighted TV
        if ( length(y) ~= length(y(:)) )
            error 'Currently only 1-dimensional signals are accepted for weighted TV'
        end
        % Currently only L1 norm is amenable for weighted TV
        if ( p ~= 1 )
            error 'Currently only L1 norm is accepted for weighted TV'
        end
        % Launch weighted solver
        x = TVL1Weighted_tautString(y,lambda);
        info = []; % Non-iterative method: no info available
        
        return;
    end

    % Choose an algorithm depending on the dimension of the input
    if isvector(y)
        % 1-dimensional solver
        [x,in] = solveTV(y, lambda, p);
        info.iters = in(1);
        info.stop = in(2);
        info.rc = in(3);
    else if length(size(y)) >= 2
        % General D-dimensional solver
        [x,in] = solveTVgen(y,lambda*ones(ndims(y),1),[1:ndims(y)],p*ones(ndims(y),1),threads, mit);
        info.iters = in(1);
        info.stop = in(2);
    end
end

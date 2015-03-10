%
% install(nopar)    proxTV installation
%
% Compiles the required mex -v files for the proxTV library.
%
% Inputs:
%   - [mode]: compilation mode
%       0: parallelization with OpenMP (default)
%       1: no parallelization
%       2: parallelization with mwlpack
%       3: debug mode
%
function install(nopar)
    % Check arguments.
    if (~exist('nopar', 'var')), nopar=0; end

    disp('Installing proxTV...');

    if(nopar == 0)
        % Compile C modules
        cd matlab
        mex -v -c -cxx CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp -fPIC' LDFLAGS='$LDFLAGS -fopenmp' ...
            ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
            ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
            ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
            ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp

        % Compile mex -v interfaces
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_condat.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_johnson.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_PN.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_tautString.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted_tautString.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2_morec.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2_PGc.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2_morec2.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVgen.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_DR.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_PD.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_CondatChambollePock.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_Yang.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2DL1W.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV3D_Yang.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVND_PDR.cpp %FIXME: add to rest of options
        mex -v -cxx -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVp_GPFW.cpp
        mex -v -cxx -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV.cpp
    elseif (nopar == 1)
        disp('WARNING: parallelization disabled');
        % Compile C modules
        mex -v -c -cxx TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex -v interfaces
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV1_condat.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV1_johnson.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV1_PN.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV1_tautString.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' TVL1Weighted.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' TVL1Weighted_tautString.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2_morec.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2_PGc.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2_morec2.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTVgen.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2D_DR.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2D_PD.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2D_CondatChambollePock.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2D_Yang.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV2DL1W.cpp
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 '*.o' solveTV3D_Yang.cpp
        mex -v -cxx -llapack -lm CXXOPTIMFLAGS=-O3 solveTVp_GPFW.cpp
        mex -v -cxx -llapack -lm CXXOPTIMFLAGS=-O3 solveTV.cpp
    elseif (nopar == 2)
        disp('Installing multisrc mwlpack version');
        % Compile C modules
        mex -v -c -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex -v interfaces
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_condat.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_johnson.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_PN.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_tautString.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted_tautString.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2_morec2.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2DL1W.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVgen.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVp_GPFW.cpp
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV.cpp
    else
        disp('WARNING: Installing in debug mode');
        % Compile C modules
        mex -v -c -DDEBUG -cxx CXXOPTIMFLAGS=-g CXXFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex -v interfaces
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_condat.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_johnson.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_PN.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV1_tautString.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' TVL1Weighted_tautString.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2_morec2.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVgen.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_DR.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2DL1W.cpp
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV2D_Condat.cpp
        mex -v -DDEBUG -cxx -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTVp_GPFW.cpp
        mex -v -DDEBUG -cxx -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '*.o' solveTV.cpp
    end
    cd ..

    % Add relevant proxTV folders to Matlab's path
    here = pwd();
    addpath([here,'/matlab']);
    addpath([here,'/demos']);
    savepath();

    disp('proxTV successfully installed.');
end

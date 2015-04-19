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
    if ~exist('nopar', 'var'), nopar=0;
    else 
        if ~isa(nopar, 'double'), error('Input argument must be numeric'); end
    end

    disp('Installing proxTV...');
    cd matlab

    if nopar == 0
        % Compile C modules
        mex -v -c -cxx CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp -fPIC' LDFLAGS='$LDFLAGS -fopenmp' ...
            ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
            ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
            ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
            ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp

        % Compile mex -v interfaces
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_PGc.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_DR.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_PD.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_CondatChambollePock.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_Yang.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV3D_Yang.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVND_PDR.cpp "*.o" %FIXME: add to rest of options
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    elseif nopar == 1
        disp('WARNING: parallelization disabled');
        % Compile C modules
        mex -v -c -cxx CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fPIC' ...
            ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
            ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
            ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
            ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp
            
        % Compile mex -v interfaces
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_condat.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_johnson.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_PN.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_tautString.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 TVL1Weighted.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 TVL1Weighted_tautString.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_morec.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_PGc.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_morec2.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTVgen.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_DR.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_PD.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_CondatChambollePock.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_Yang.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2DL1W.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV3D_Yang.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTVND_PDR.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTVp_GPFW.cpp "*.o"
        mex -v -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV.cpp "*.o"
    elseif nopar == 2
        disp('Installing multisrc mwlpack version');
        % Compile C modules
        mex -v -c -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS='$CXXFLAGS -fopenmp -fPIC' LDFLAGS='$LDFLAGS -fopenmp' ...
            ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
            ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
            ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
            ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp

        % Compile mex -v interfaces
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_PGc.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_DR.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_PD.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_CondatChambollePock.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_Yang.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV3D_Yang.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVND_PDR.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -v -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    else
        disp('WARNING: Installing in debug mode');
        % Compile C modules
        mex -v -c -DDEBUG -cxx CXXOPTIMFLAGS=-g CXXFLAGS='$CXXFLAGS -fopenmp -fPIC' LDFLAGS='$LDFLAGS -fopenmp' ...
            ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
            ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
            ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
            ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp
        
        % Compile mex -v interfaces
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_PGc.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_DR.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_PD.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_CondatChambollePock.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_Yang.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV3D_Yang.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVND_PDR.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -v -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    end
    cd ..

    % Add relevant proxTV folders to Matlab's path
    here = pwd();
    addpath([here,'/matlab']);
    addpath([here,'/matlab/demos']);
    savepath();

    disp('proxTV successfully installed.');
end

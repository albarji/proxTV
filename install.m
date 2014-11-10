%
% install(nopar)    proxTV installation
%
% Compiles the required mex files for the proxTV library.
%
% Inputs:
%   - [mode]: compilation mode
%       0: parallelization with OpenMP (default)
%       1: no parallelization
%       2: parallelization with mwlpack
%       3: debug mode
%
function install(nopar)
    % Check arguments 
    if (~exist('nopar', 'var')), nopar=0; end
    
    disp('Installing proxTV...');
    cd src
    
    if(nopar == 0)
        % Compile C modules
        mex -c -cxx CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex interfaces
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_PGc.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_DR.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_PD.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_CondatChambollePock.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_Yang.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV3D_Yang.cpp "*.o"
        mex -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVND_PDR.cpp "*.o" %FIXME: add to rest of options
        mex -cxx  -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -cxx  -llapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    elseif (nopar == 1)
        disp('WARNING: parallelization disabled'); 
        % Compile C modules
        mex -c -cxx TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex interfaces
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_condat.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_johnson.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_PN.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV1_tautString.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 TVL1Weighted.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 TVL1Weighted_tautString.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_morec.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_PGc.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2_morec2.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTVgen.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_DR.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_PD.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_CondatChambollePock.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2D_Yang.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV2DL1W.cpp "*.o"
        mex -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-O3 solveTV3D_Yang.cpp "*.o"
        mex -cxx -llapack -lm CXXOPTIMFLAGS=-O3 solveTVp_GPFW.cpp "*.o"
        mex -cxx -llapack -lm CXXOPTIMFLAGS=-O3 solveTV.cpp "*.o"
    elseif (nopar == 2)
        disp('Installing multicore mwlpack version');
        % Compile C modules
        mex -c -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex interfaces
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -cxx -lmwlapack -lm CXXOPTIMFLAGS=-O3 CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    else
        disp('WARNING: Installing in debug mode');
        % Compile C modules
        mex -c -DDEBUG -cxx CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVgenopt.cpp TVL1opt.cpp TVL1Wopt.cpp TVL2opt.cpp TVLPopt.cpp TV2Dopt.cpp TV2DWopt.cpp TVNDopt.cpp LPopt.cpp utils.cpp condat_fast_tv.cpp johnsonRyanTV.cpp
        % Compile mex interfaces
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_condat.cpp "*.o"
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_johnson.cpp "*.o"
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_PN.cpp "*.o"
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV1_tautString.cpp "*.o"
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted.cpp "*.o"
        mex -DDEBUG -cxx -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TVL1Weighted_tautString.cpp "*.o"
        mex -DDEBUG -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2_morec2.cpp "*.o"
        mex -DDEBUG -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVgen.cpp "*.o"
        mex -DDEBUG -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_DR.cpp "*.o"
        mex -DDEBUG -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2DL1W.cpp "*.o"
        mex -DDEBUG -cxx  -lblas -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV2D_Condat.cpp "*.o"
        mex -DDEBUG -cxx  -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTVp_GPFW.cpp "*.o"
        mex -DDEBUG -cxx  -llapack -lm CXXOPTIMFLAGS=-g CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" solveTV.cpp "*.o"
    end
    cd ..
    
    % Add relevant proxTV folders to Matlab's path
    here = pwd();
    addpath([here,'/src']);
    addpath([here,'/demos']);
    savepath();
    
    disp('proxTV successfully installed.');
end

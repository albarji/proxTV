%
% install_octave(nopar)    proxTV installation
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
function install_octave(nopar)
    % Check arguments.
    if ~exist('nopar', 'var'), nopar=0;
    else 
        if ~isa(nopar, 'double'), error('Input argument must be numeric'); end
    end

    disp('Installing proxTV...');
    cd matlab

    %% compilation flags 
    [~, CXXFLAGS] = system('mkoctfile -p CXXFLAGS');
    [~, LDFLAGS] = system('mkoctfile -p LDFLAGS');
    % some versions introduces a newline character (10)
    % in the output of 'system'; this must be removed
    if CXXFLAGS(end)==10, CXXFLAGS = CXXFLAGS(1:end-1); end
    if LDFLAGS(end)==10, LDFLAGS = LDFLAGS(1:end-1); end
    CXXFLAGSorig = CXXFLAGS;
    LDFLAGSorig = LDFLAGS;

    CXXFLAGS = [CXXFLAGS ' -fPIC -DOCTAVE'];

    if nopar == 0
        CXXFLAGS = [CXXFLAGS ' -O3 -fopenmp'];
        LDFLAGS = [LDFLAGS ',-fopenmp'];
        setenv('CXXFLAGS', CXXFLAGS);
        setenv('LDFLAGS', LDFLAGS);
    elseif nopar == 1
        disp('WARNING: parallelization disabled');
        CXXFLAGS = [CXXFLAGS ' -O3'];
        setenv('CXXFLAGS', CXXFLAGS);
    elseif nopar == 2
        disp('Installing multisrc mwlpack version');
        CXXFLAGS = [CXXFLAGS ' -O3 -fopenmp'];
        LDFLAGS = [LDFLAGS ',-fopenmp'];
    else
        disp('WARNING: Installing in debug mode');
        CXXFLAGS = [CXXFLAGS ' -g -fopenmp'];
        LDFLAGS = [LDFLAGS ',-fopenmp'];
    end

    setenv('CXXFLAGS', CXXFLAGS);
    setenv('LDFLAGS', LDFLAGS);

    %% Compile C modules
    mex -v -c ../src/TVgenopt.cpp ../src/TVL1opt.cpp ../src/TVL1Wopt.cpp ...
        ../src/TVL2opt.cpp ../src/TVLPopt.cpp ../src/TV2Dopt.cpp ...
        ../src/TV2DWopt.cpp ../src/TVNDopt.cpp ../src/LPopt.cpp ...
        ../src/utils.cpp ../src/condat_fast_tv.cpp ../src/johnsonRyanTV.cpp ...
        ../src/TVL1opt_tautstring.cpp ../src/TVL1opt_hybridtautstring.cpp ...
        ../src/TVL1opt_kolmogorov.cpp

    fprintf('\n');

    %% Compile MEX
    solvers = {'solveTV1_condat', 'solveTV1_condattautstring', ...
        'solveTV1_johnson', 'solveTV1_PN', 'solveTV1_linearizedTautString', ...
        'solveTV1_classicTautString', 'solveTV1_hybridTautString', ...
        'solveTV1_kolmogorov', 'TVL1Weighted', 'TVL1Weighted_tautString', ...
        'solveTV2_morec', 'solveTV2_PGc', 'solveTV2_morec2', 'solveTVgen', ...
        'solveTV2D_DR', 'solveTV2D_PD', 'solveTV2D_CondatChambollePock', ...
        'solveTV2D_Yang', 'solveTV2D_Kolmogorov', 'solveTV2DL1W', ...
        'solveTV3D_Yang', 'solveTVND_PDR', 'solveTVp_GPFW', 'solveTV'};
    for solver=solvers
        if nopar == 0 || nopar == 1
            eval(['mex -v -lblas -llapack -lm ' solver{:} '.cpp "*.o"']);
        elseif nopar == 2
            eval(['mex -v -lmwlapack -lm ' solver{:} '.cpp "*.o"']);
        else
            eval(['mex -v -DDEBUG -lblas -llapack -lm' solver{:} ...
                '.cpp "*.o"']);
        end
        system(['rm -f ' solver{:} '.o']); % some versions keep the object file
    end
    clear *.mex
    system('rm -f *.o');

    cd ..

    % Add relevant proxTV folders to Matlab's path
    % here = pwd();
    % addpath([here,'/matlab']);
    % addpath([here,'/matlab/demos']);
    % savepath();

    setenv('CXXFLAGS', CXXFLAGSorig);
    setenv('LDFLAGS', LDFLAGSorig);

    disp('proxTV successfully installed.');
end

addpath(genpath(pwd));

% --- paths ---
root = pwd;
addpath(genpath(fullfile(root,'data')));
addpath(genpath(fullfile(root,'src'))); 

% --- compile MEX only if missing ---
mexfile = fullfile(root,'src','solver_functions', ['ks_nonlinear_mex.' mexext]);

% --- build if needed ---
if exist(mexfile,'file') ~= 3
    fprintf('Compiling C++ codes...\n');
    mex -O ...
        CXXFLAGS="\$CXXFLAGS -fopenmp" ...
        LDFLAGS="\$LDFLAGS -fopenmp" ...
        -outdir src/solver_functions ...
        src/solver_functions/ks_nonlinear_mex.cpp
end

% --- run solver ---
main_2DKS(1e-5,256,5,5,1,2.11,2.11,.02,-3.0,-3.0,1,'optimize','IC',1e-6,0.0);
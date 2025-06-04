%function main_2DKS_v2_0
tic

%%% choose switches %%%
run = 'IC'; % switch to 0, 'L', 'N', 'dt', 'T', 'IC', 'kappa'
test_parameter = 0; % do not alter
diagnostics = 1; % evaluating dynamics of Fourier modes
method  = 'imexrk4vec2a'; % time-stepping scheme
%IC = 'mn1'; % initial condition

%%% choose parameter ranges %%%
timewindow = linspace(0,60,7); % for time-stepping analysis
timewindow(1) = 1;
%L_scale = linspace(11,19,9)/10; % for dynamical behavior analysis
%L_scale = [ 1.1 , 1.4 , 1.5 , 1.9 , 3.2 , 5.2 , 10.2]; % for dynamical behavior analysis
L_scale = 1.00:0.01:2.00;
timestep = 10.^(-linspace(1,4,7)); % for temporal convergence analysis
gridsize = 10*4*linspace(3,21,7); % for spatial convergence analysis
initialcondition = { 'sinL' };

%%% choose default parameters %%%
%{
L_s1 = L_scale(6); % length-scale parameter in dim 1
L_s2 = L_scale(4); % length-scale parameter in dim 2
dt = timestep(5); % length of time-step
T = timewindow(2); % time window
N = gridsize(3); % number of grid points 
save_each = 100; % number of iterations between saved timepoints
%}

for lscale = 1:1%length(L_scale)

    %%% choose temporary parameters %%%
    %
    L_s1 = L_scale(lscale); % length-scale parameter in dim 1
    L_s2 = L_s1; % length-scale parameter in dim 2
    dt = 1e-2; % length of time-step
    T = 1000; % time window
    N = 32; % number of grid points 
    save_each = 1/dt; % number of iterations between saved timepoints - use T*10 to save 100, 1/dt to save 1 T
    %}

    switch run % set which test to run
        case 'L'
            test_parameter = L_scale; % (dynamical behavior)
        case 'N'
            test_parameter = gridsize; % (spatial convergence)
        case 'dt'
            test_parameter = timestep; % (temporal convergence)
        case 'T'
            test_parameter = timewindow; % (temporal convergence)
        case 'IC'
            test_parameter = initialcondition; % 
        case 'kappa'
            test_parameter = 1; % (dummy variable)
            save_each = 1; % save all timesteps for adjoint solver
    end

    %
    for k = 1 : length(test_parameter) % length('x') indicates 'x' testing

        method  = 'imexrk4vec2a'; % time-stepping scheme
    
        %%% set variable parameters %%%
        switch run 
            case 'L'
                L_s2 = L_scale(k); % length-scale parameter in dim 2 
            case 'N'
                N = gridsize(k); % number of grid points 
            case 'dt'
                dt = timestep(k); % length of time-step
            case 'T'
                T = timewindow(k); % length of simulation time window
            case 'IC'
                IC = strjoin(initialcondition(k),'');
        end
        
        %%% solve PDE problem in time %%%
        %
        tic
        [ v_n , u_n ] = DirectSolve_2DKS_v2_0(IC,method,N,L_s1,L_s2,dt,T,save_each);
        time = toc;
        toc
        %}

        %%% save/inspect solution %%%
        switch run 
            case {'L','N','dt','T','IC'}                
                %
                method = IC;
                save_2DKSsolution('time_evolution', u_n, v_n, time, method, dt, T, N, L_s1, L_s2); % save solution
                %{
                [u_n, v_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2); % load solution
                %}
                %plot2DKS(v_n , u_n, 'initial', method, N, dt, T, L_s1, L_s2); % save/inspect initial state
                %plot2DKS(v_n , u_n, 'terminal', method, N, dt, T, L_s1, L_s2); % save/inspect terminal state
                plot2DKS(v_n , u_n, 'diagnostics', method, N, dt, T, L_s1, L_s2); % save/inspect dynamical characteristics
                %plot2DKS(v_n , u_n, 'gif', method, N, dt, T, L_s1, L_s2); % save/inspect surface time evolution                
                close all
                %
        end

    end
    %

    %{
    switch run 
        case 'N' % spatial convergence: error analysis and computational time
            [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, method, dt, T, L_s1, L_s2);
            save_measures('spatial', error_2, error_inf, comptime, method, 0, dt, T, L_s1, L_s2);
        case 'dt' % temporal convergence: error analysis and computational time
            [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, method, N, T, L_s1, L_s2);
            save_measures('temporal', error_2, error_inf, comptime, method, N, 0, T, L_s1, L_s2);
    end
    %}

end

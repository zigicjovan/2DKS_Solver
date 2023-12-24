function main_2DKS_v2_0

    %%% choose switches %%%
    run = 'L'; % switch to 0, 'L', 'N', 'dt', 'T'
    test_parameter = 0; % do not alter
    method  = 'imexrk4vec2a'; % time-stepping scheme
    IC = 'sin'; % initial condition

    %%% choose parameter ranges %%%
    timewindow = linspace(0,60,7); % for time-stepping analysis
    timewindow(1) = 1;
    L_scale = linspace(0,20,11)/10; % for dynamical behavior analysis
    L_scale(1) = 10^(-1.5);
    timestep = 10.^(-linspace(1,4,7)); % for temporal convergence analysis
    gridsize = 10*4*linspace(3,21,7); % for spatial convergence analysis

    %%% choose default parameters %%%
    L_s1 = 1.8; % length-scale parameter in dim 1
    L_s2 = L_scale(4); % length-scale parameter in dim 2
    dt = timestep(5); % length of time-step
    T = timewindow(2); % time window
    N = gridsize(3); % number of grid points 
    save_each = 100; % number of iterations between saved timepoints

    switch run % set which test to run
        case 'L'
            test_parameter = L_scale; % (dynamical behavior)
        case 'N'
            test_parameter = gridsize; % (spatial convergence)
        case 'dt'
            test_parameter = timestep; % (temporal convergence)
        case 'T'
            test_parameter = timewindow; % (temporal convergence)
    end

    %
    for k = 1 : length(test_parameter) % length('x') indicates 'x' testing

        close all;
    
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
        end
        
        %%% solve PDE problem in time %%%
        %
        tic
        [ v_n , u_n ] = DirectSolve_2DKS_v1_0(IC,method,N,L_s1,L_s2,dt,T,save_each);
        time = toc;
        %}
        
        %%% save/inspect solution %%%
        %
        save_2DKSsolution('time_evolution', u_n, v_n, time, method, dt, T, N, L_s1, L_s2); % save solution
        %{
        [u_n, v_n, ~] = load_2DKSsolution('time_evolution', method, dt, T, N, L_s1, L_s2); % load solution
        %}
        plot2DKS(v_n , u_n, 'initial', method, N, dt, T, L_s1, L_s2); % save/inspect initial state
        plot2DKS(v_n , u_n, 'terminal', method, N, dt, T, L_s1, L_s2); % save/inspect terminal state
        plot2DKS(v_n , u_n, 'gif_contour', method, N, dt, T, L_s1, L_s2); % save/inspect contour time evolution
        plot2DKS(v_n , u_n, 'gif', method, N, dt, T, L_s1, L_s2); % save/inspect surface time evolution
        %
    
    end
    %

    %
    switch run 
        case 'N' % spatial convergence: error analysis and computational time
            [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, method, dt, T, L_s1, L_s2);
            save_measures('spatial', error_2, error_inf, comptime, method, 0, dt, T, L_s1, L_s2);
        case 'dt' % temporal convergence: error analysis and computational time
            [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, method, N, T, L_s1, L_s2);
            save_measures('temporal', error_2, error_inf, comptime, method, N, 0, T, L_s1, L_s2);
    end
    %

end
%function main_2DKS_v2_0
tic

%%% choose switches %%%
run = 'IC';                                  % switch to 0, 'L', 'N', 'dt', 'T', 'IC', 'kappa'
test_parameter = 0;                             % do not alter
diagnostics = 1;                                % evaluating dynamics of Fourier modes
IC = 'sinL';                                    % initial condition

%%% choose parameter ranges %%%
timewindow = 2000;%linspace(20,200,19);               % for time-stepping analysis
L_scale = 1:0.01:2.8;                         % for dynamical behavior analysis
timestep = 10.^(-linspace(1,4,7));              % for temporal convergence analysis
gridsize = 10*4*linspace(3,21,7);               % for spatial convergence analysis
initialcondition = { 'sinL' };

%%% choose default parameters %%%
%{
L_s1 = L_scale(6);                              % length-scale parameter in dim 1
L_s2 = L_scale(4);                              % length-scale parameter in dim 2
dt = timestep(5);                               % length of time-step
T = timewindow(2);                              % time window
N = gridsize(3);                                % number of grid points 
save_each = 100;                                % number of iterations between saved timepoints
%}

%kappalist_48_1e2_sinL10_p30 = NaN(15,length(L_scale)*length(timewindow)); % temporary for kappa tests
%rieszlist_48_1e2_sinL_p1_T200 = NaN(length(L_scale)*length(timewindow),1); % temporary for kappa tests

for lscale = 1:length(L_scale)

    %%% choose temporary parameters %%%
    %
    L_s1 = L_scale(lscale);                     % length-scale parameter in dim 1
    L_s2 = L_s1;                                % length-scale parameter in dim 2
    dt = 1e-2;                                  % length of time-step
    T = 2000;                                    % time window
    N = 48;                                     % number of grid points 
    save_each = 2/dt;                           % number of iterations between saved timepoints - use T*10 to save 100, 1/dt to save 1 T
    %}

    switch run % set which test to run
        case 'L'
            test_parameter = L_scale;           % (dynamical behavior)
        case 'N'
            test_parameter = gridsize;          % (spatial convergence)
        case 'dt'
            test_parameter = timestep;          % (temporal convergence)
        case 'T'
            test_parameter = timewindow;        % (temporal convergence)
        case 'IC'
            test_parameter = initialcondition;  
        case 'kappa'
            test_parameter = timewindow;        % (dummy variable)
            save_each = 1;                      % save all timesteps for backward solver
            J_pert = NaN(15,1);                 % store perturbed objective functionals  
            gateaux_deriv = NaN(15,1);          % store kappa test numerators 
            kappa = NaN(15,1);                  % store kappa test results 
    end

    for k = 1 : length(test_parameter)          % length('x') indicates 'x' testing

        %%% set variable parameters %%%
        switch run 
            case 'L'
                L_s2 = L_scale(k);              % length-scale parameter in dim 2 
            case 'N'
                N = gridsize(k);                % number of grid points 
            case 'dt'
                dt = timestep(k);               % length of time-step
            case 'T'
                T = timewindow(k);              % length of simulation time window
            case 'IC'
                IC = strjoin(initialcondition(k),'');
            case 'kappa'
                T = timewindow(k);              % length of simulation time window
                save_each = 1;                      % save all timesteps for backward solver
        end
        
        %%% solve PDE problem in time %%%
        disp(['Solving forward equation for L = ' num2str(L_s1) ', T = ' num2str(T)])
        tic
        [ v_n , u_n ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,0);
        time = toc;
        toc

        %%% save/inspect solution %%%
        switch run 
            case {'L','N','dt','T','IC'}                
                save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, L_s1, L_s2);               % save solution
                %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2);              % load solution
                %plot2DKS(u_n, 'initial', IC, N, dt, T, L_s1, L_s2, 0);                                 % save/inspect initial state
                %plot2DKS(u_n, 'terminal', IC, N, dt, T, L_s1, L_s2, 0);                                % save/inspect terminal state
                %plot2DKS(u_n, 'diagnostics', IC, N, dt, T, L_s1, L_s2, 0);                             % save/inspect dynamical characteristics
                plot2DKS(u_n, 'norms', IC, N, dt, T, L_s1, L_s2, 0);                                   % save/inspect dynamical characteristics
                %plot2DKS(u_n, 'gif', IC, N, dt, T, L_s1, L_s2, 0);                                     % save/inspect time evolution 
                toc
            case 'kappa' 
                L1 = 2*pi*L_s1;
                L2 = 2*pi*L_s2;
                J_init = sum( u_n(:,end) .* conj(u_n(:,end)) )*(L1*L2)/N^2;                             % initial objective functional (L^2 inner product of terminal forward state)        
                disp('Solving backward equation...')
                [v_adj, u_adj, u_pertIC] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,v_n);              % solve adjoint equation
                toc
                gat_riesz = sum( u_adj(:,1) .* conj(u_pertIC) )*(L1*L2)/N^2;                            % kappa test denominator 
                rieszlist_48_1e2_sinL_p1_T200((lscale-1)*length(timewindow)+k,1) = gat_riesz;                 % temporary for kappa tests
                gradJ = v_adj(:,1);                                                                     % objective gradient
                save_each = 1/dt;                                                                       % release memory
                clear v_n v_adj u_adj                                                                   % release memory
                disp('Solving perturbed equations...')
                for i = -15:1:-1
                    eps = 10^(i);                                                                       % define perturbation 
                    [ ~ , u_pert ] = solve_2DKS(IC,'kappa',N,L_s1,L_s2,dt,T,save_each,eps);             % solve perturbed forward equation
                    J_pert(i+16,1) = sum( u_pert(:,end) .* conj(u_pert(:,end)) )*(L1*L2)/N^2;           % perturbed objective functional
                    gateaux_deriv(i+16,1) = (J_pert(i+16,1) - J_init)/eps;                              % kappa test numerator 
                    kappa(i+16,1) = gateaux_deriv(i+16,1)/gat_riesz;                                    % kappa test numerator 
                    toc
                end
                clear v_pert u_pert                                                                     % release memory
                plot2DKS(u_n, 'kappa', IC, N, dt, T, L_s1, L_s2, kappa);                                % save/inspect kappa test figure
                %kappalist_48_1e2_sinL10_p30(:,(lscale-1)*length(timewindow)+k) = kappa;
        end

        close all                                                                                       % close any open figures

    end

    switch run 
        case 'N'                                                                                        % spatial convergence: error analysis and computational time
            [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, IC, dt, T, L_s1, L_s2);
            save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, L_s1, L_s2);
        case 'dt'                                                                                       % temporal convergence: error analysis and computational time
            [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, IC, N, T, L_s1, L_s2);
            save_measures('temporal', error_2, error_inf, comptime, meICthod, N, 0, T, L_s1, L_s2);
    end

end

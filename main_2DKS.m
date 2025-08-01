%function main_2DKS
tic

%%% choose test %%%
run = 'optimize';                               % switch to 'L', 'N', 'dt', 'IC', 'kappa', 'energygrowth', 'optimize'
Ntime_save_max = 10000;                         % choose maximum number of samples per data file

%%% choose parameter testing ranges %%%
L_scale = 2.36;%[1.1,1.4,1.5,2.2,3.2,5.2,10.2];                                % domain sizes
timestep = .005;                                % time-step sizes
gridsize = 48;                                  % grid sizes
timewindow = 5.5;                                % time windows
initialcondition = {'s1'};                      % initial conditions
kappapert = {'s1'};                             % perturbation functions
L_target = 2.36;                                % domain sizes of interest

%%% choose default parameters %%%
L_s1 = L_scale(1);                              % length-scale parameter in dim 1
L_s2 = L_scale(1);                              % length-scale parameter in dim 2
dt = timestep(1);                               % length of time-step
N = gridsize(1);                                % number of grid points
T = timewindow(1);                              % length of simulation time window
IC = strjoin(initialcondition(1),'');           % initial condition
save_each = 1;                                  % number of iterations between saved timepoints - 1/dt to save each 1 T

numberoftests = length(initialcondition)*length(kappapert)*length(timewindow)*length(L_scale);
testcounter = 0;

for init = 1 : length(kappapert)

    pertIC = strjoin(kappapert(init),'');

    for choros = 1 : length(L_scale)
    
        L_s1 = L_scale(choros);                         % adjust length-scale parameter in dim 1
        L_s2 = L_s1;                                    % adjust length-scale parameter in dim 2
    
        for chronos = 1 : length(timewindow)
    
            T = timewindow(chronos);                    % adjust simulation time window
    
            switch run                                  % set which test to run
                case 'L'
                    test_parameter = L_scale;           % (dynamical behavior)
                case 'N'
                    test_parameter = gridsize;          % (spatial convergence)
                case 'dt'
                    test_parameter = timestep;          % (temporal convergence)
                case {'IC', 'energygrowth','kappa','optimize'}
                    test_parameter = initialcondition;  
            end
        
            for k = 1 : length(test_parameter)          % length('x') indicates 'x' testing
        
                %%% set variable parameters %%%
                switch run 
                    case 'L'
                        L_s2 = L_scale(k);                      % length-scale parameter in dim 2 
                    case 'N'
                        N = gridsize(k);                        % number of grid points 
                    case 'dt'
                        dt = timestep(k);                       % length of time-step
                    case {'IC', 'energygrowth'}
                        IC = strjoin(initialcondition(k),'');   % choice of initial condition
                    case {'kappa','optimize'}
                        save_each = 1;
                        IC = strjoin(initialcondition(k),'');   % choice of initial condition
                end
                
                testcounter = testcounter + 1;
                disp(['Test ' num2str(testcounter) ' of ' num2str(numberoftests) ])
                %%% solve forward-time PDE problem %%%
                switch run 
                    case {'L','N','dt','IC','kappa','optimize'} 
                        disp(['Solving forward-time problem for L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt)])
                        tic
                        [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                        time = toc;
                        toc
                end
        
                %%% save/inspect solution %%%
                switch run 
                    case 'energygrowth'
                        if exist('l2norms_avg','var') == 0
                            l2norms_avg = NaN(length(L_scale),length(test_parameter)+1);
                            l2norms_mode = NaN(length(L_scale),length(test_parameter)+1);
                            l2norms_guess = NaN(T/dt,numberoftests);
                            l2norms_opt = NaN(T/dt,numberoftests);
                        end
                        %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2, 0);                   % load solution
                        [u_normL2, ~] = load_2DKSsolution('normL2', IC, dt, T, N, L_s1, L_s2, 0, 0);                       % load solution
                        l2norms_guess(:,testcounter) = u_normL2;
                        %[u_normL2, ~] = load_2DKSsolution('normL2', 'optimized', dt, T, N, L_s1, L_s2, 0, IC);                       % load solution
                        l2norms_opt(:,testcounter) = u_normL2;
                        u_normL2rd = round(u_normL2,1);
                        l2norms_mode(choros,k) = mode(u_normL2rd);
                        l2norms_avg(choros,k) = mean(u_normL2(ceil(end/2):end,1));
                        if k == length(test_parameter) 
                            l2norms_avg(choros,k+1) = std(l2norms_avg(choros,1:k));
                            l2norms_mode(choros,k+1) = std(l2norms_mode(choros,1:k));
                        end
                    case {'L','N','dt','IC'}                
                        %save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, L_s1, L_s2,0,0);                      % save solution
                        %plot_2DKS(save_each, 'initial', IC, N, dt, T, L_s1, L_s2,Ntime_save_max, 0,0);                                % save/inspect initial state
                        %plot_2DKS(save_each, 'terminal', IC, N, dt, T, L_s1, L_s2,Ntime_save_max, 0,0);                               % save/inspect terminal state
                        %plot_2DKS(save_each, 'diagnostics', IC, N, dt, T, L_s1, L_s2,Ntime_save_max, 0,0);                            % save/inspect dynamical characteristics
                        plot_2DKS(save_each, 'norms', IC, N, dt, T, L_s1, L_s2,Ntime_save_max, 0,0);                                  % save/inspect dynamical characteristics
                        %plot_2DKS(save_each, 'gif', IC, N, dt, T, L_s1, L_s2,Ntime_save_max, 0,100);                                     % save/inspect time evolution 
                        toc
                        close all                                                                                       % close any open figures
                    case 'kappa' 
                        %pertIC = IC;
                        if testcounter > 1
                            [kappa,gat_riesz,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,kappalist,rieszlist,Ntime_save_max);
                        else
                            [kappa,gat_riesz,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,0,0,Ntime_save_max);
                        end               
                        delete_2DKSsolution('forward', IC, dt, T, N, L_s1, L_s2, Ntime_save_max,0);
                        delete_2DKSsolution('backward', IC, dt, T, N, L_s1, L_s2, Ntime_save_max,0);
                    case 'optimize' 
                        [J_opt, J_history , u_TC_opt , u_IC_opt] = optimize_2DKS(IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,u_IC,Ntime_save_max);
                        %plot_2DKS(save_each, 'gif', 'optimized', N, dt, T, L_s1, L_s2,Ntime_save_max,IC,0);   
                        plot_2DKS(save_each, 'diagnostics', 'optimized', N, dt, T, L_s1, L_s2,Ntime_save_max,IC ,0);
                        plot_2DKS(save_each, 'initial', IC, N, dt, T, L_s1, L_s2,Ntime_save_max,IC,0);
                        plot_2DKS(save_each, 'terminal', IC, N, dt, T, L_s1, L_s2,Ntime_save_max,IC,0);
                        plot_2DKS(save_each, 'initial', 'optimized', N, dt, T, L_s1, L_s2,Ntime_save_max,IC,0);
                        plot_2DKS(save_each, 'terminal', 'optimized', N, dt, T, L_s1, L_s2,Ntime_save_max,IC,0);
                        disp(['Solved optimization problem for L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt)])
                        disp(['Initial objective functional value: ' num2str(J_history(1,1))])
                        disp(['Optimal objective functional value: ' num2str(J_history(end,1))])
                        disp(['Number of iterations: ' num2str(length(J_history)-1)])
                        save_2DKSsolution('optimal', u_IC_opt, 0, 0, IC, dt, T, N, L_s1, L_s2, 1, 0); % save solution to machine
                        close all
                end
            end
        end
        switch run 
            case 'kappa'    % designed for 5 or 10 tests only
                %plot_measures('kappa', dt, pertIC, N, timewindow, L_s1, testcounter, length(timewindow));
        end
    end
end

switch run 
    case 'energygrowthXX'     % designed for comparing two initial conditions only 
        plot_measures('energygrowth', L_scale, initialcondition, l2norms_mode, T, L_target, l2norms_avg);
        save_measures('energygrowth', l2norms_mode, l2norms_avg, 0, initialcondition, 0, L_scale, T, L_target, 0);
    case 'N'                % spatial convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('spatial', dt, IC, gridsize, T, L_s1, L_s2);
        save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, L_s1, L_s2);
    case 'dt'               % temporal convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('temporal', timestep, IC, N, T, L_s1, L_s2);
        save_measures('temporal', error_2, error_inf, comptime, IC, N, 0, T, L_s1, L_s2);
end
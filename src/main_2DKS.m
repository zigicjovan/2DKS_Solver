function main_2DKS(dtc,Nc,Kstart,Kend,Knum,ell1start,ell1end,ell1gap,ell2start,ell2end,ell2gap,Tstart,Tend,Tnum,runc,continuationc,guessc,tolc,optTc)

% --- paths ---
root = pwd;
addpath(genpath(fullfile(root,'data')));
addpath(genpath(fullfile(root,'src'))); 

% --- compile MEX only if missing ---
ensure_mex('src/solver_functions/ks_nonlinear_fwd.cpp','src/solver_functions');
ensure_mex('src/solver_functions/ks_nonlinear_bwd.cpp','src/solver_functions');

tic

%%% choose test settings %%%
restart = 1;                                    % binary switch to generate new test counters
optfigs = 1;                                    % generate optimization diagnostic figures 
Ntime_save_max = 10;                            % choose maximum number of samples per data file
run = runc;                               	    % switch to 'optimize', 'plotOptIC', 'energygrowth', 'L', 'N', 'dt', 'IC', 'kappa'
continuation = continuationc;                   % 'IC' for optimal IC from file, 'off' to generate new data
optmethod = 'RCG';                              % RCG, RG, or RCGd5 (start after 5th iter)
timestep = dtc;                                 % time-step sizes
gridsize = Nc;                                  % grid sizes
tol = tolc;                                     % set optimization tolerance critera
optT = 10^(optTc);                              % T parameter of optimal IC for asymptotic simulations

%%% choose parameter testing ranges %%%
initialKmagnitude = logspace(Kstart,Kend,Knum);                        % initial L^2 energy magnitudes
L1_scale = ell1start:ell1gap:ell1end;                        % domain sizes
L2_scale = ell2start:ell2gap:ell2end;
timewindow = logspace(Tstart,Tend,Tnum);         % time windows
initialcondition = guessc; % initial conditions

%timewindow = timewindow(2);
switch run
    case 'kappa'
        kappapert = {'stg30'};                          % perturbation functions
        L_target = [2.36,2.78];                         % domain sizes of interest
end

%%% choose default parameters %%%
K = initialKmagnitude(1);                       % magnitude of initial condition L^2 energy, use K=0 for default IC norm
L_s1 = L1_scale(1);                              % length-scale parameter in dim 1
L_s2 = L2_scale(1);                              % length-scale parameter in dim 2
T = timewindow(1);                              % length of simulation time window
IC = strjoin(initialcondition(1),'');           % initial condition
dt = timestep(1);                               % length of time-step
N = gridsize(1);                                % number of grid points
save_each = 1;                                  % number of iterations between saved timepoints - 1/dt to save each 1 T

switch optmethod
    case 'RCG'
        RCGon = 1;                              % RCG switch for media files
    otherwise
        RCGon = 0;
end

if restart == 1
    numberoftests = length(initialcondition)*length(initialKmagnitude)*length(timewindow)*length(L1_scale)*length(L2_scale);
    testcounter = 0;
    testrow = 0;
    
    Joptdata = NaN(length(initialKmagnitude)*length(timewindow)*length(L1_scale)*length(L2_scale)*length(initialcondition),3+3);
    Jinitdata = NaN(length(initialKmagnitude)*length(timewindow)*length(L1_scale)*length(L2_scale)*length(initialcondition),3+1);
end

for energy_i = 1 : length(initialKmagnitude)

    switch run
        case 'kappa'
            pertIC = strjoin(kappapert(energy_i),'');
    end
    K = initialKmagnitude(energy_i);                 % adjust initial L^2 energy magnitude

    for domain_i = 1 : length(L1_scale)
    
        L_s1 = L1_scale(domain_i);                         % adjust length-scale parameter in dim 1
        L_s2 = L2_scale(1);                                    % adjust length-scale parameter in dim 2
    
        %% start window block
        for window_i = 1 : length(timewindow)
    
            T = timewindow(window_i);                    % adjust simulation time window
    
            switch run                                  % set which test to run
                case 'L'
                    test_parameter = L1_scale;           % (dynamical behavior)
                case 'N'
                    test_parameter = gridsize;          % (spatial convergence)
                case 'dt'
                    test_parameter = timestep;          % (temporal convergence)
                case {'IC', 'energygrowth','kappa','optimize','plotOptIC'}
                    test_parameter = initialcondition;  
            end
        
            testrow = testrow + 1;

            for param_i = 1 : length(test_parameter)          % length('x') indicates 'x' testing
        
                %%% set variable parameters %%%
                switch run 
                    case 'L'
                        L_s2 = L1_scale(param_i);                      % length-scale parameter in dim 2 
                    case 'N'
                        N = gridsize(param_i);                        % number of grid points 
                    case 'dt'
                        dt = timestep(param_i);                       % length of time-step
                    case {'IC', 'energygrowth','optimize','plotOptIC'}
                        IC = strjoin(initialcondition(param_i),'');   % choice of initial condition
                    case {'kappa'}
                        save_each = 1;
                        IC = strjoin(initialcondition(param_i),'');   % choice of initial condition
                end
                
                testcounter = testcounter + 1;
                disp(['Test ' num2str(testcounter) ' of ' num2str(numberoftests) ])
                %%% solve forward-time PDE problem %%%
                switch run 
                    case {'plotOptIC'}
                        initIC = IC;
                        optdt = dt;
                        %optT = optT;
                        optN = N;
                        optK = K;
                        optL_s1 = L_s1;
                        optL_s2 = L_s2;
                        opttol = tol; 
                        [ u_IC , ~ ] = load_2DKSsolution('optimal', initIC, optdt, optT, optN, optK, optL_s1, optL_s2, [opttol optT], 0);
                        tic
                        disp(['Using chosen optimal IC for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                        newopt = 0;
                        plot_2DKS(save_each, 'plotOptIC', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,u_IC,u_IC,[tol,RCGon,100,newopt]);
                        fprintf('Simulation run complete.\n')
                    case {'L','N','dt','IC','kappa'} 
                        switch continuation
                            case 'IC' % choose which parameters define optimized initial data
                                initIC = IC;
                                optdt = dt;
                                %optT = optT;
                                optN = N;
                                optK = K;
                                optL_s1 = L_s1;
                                optL_s2 = L_s2;
                                opttol = tol; 
                                [ u_IC , ~ ] = load_2DKSsolution('optimal', initIC, optdt, optT, optN, optK, optL_s1, optL_s2, [opttol optT], 0);
                                tic
                                disp(['Using chosen optimal IC for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                %solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                originalIC = IC;
                                IC = 'optimized';
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC,originalIC);
                                disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                            case 'off'
                                disp(['Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                tic
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                        end
                    case 'optimize'
                        originalIC = IC;
                        switch continuation
                            case 'IC'
                                try
                                    [ u_IC , ~ ] = load_2DKSsolution('optimal', IC, dt, T, N, K, L_s1, L_s2, [tol T], 0);
                                    tic
                                    disp(['Continuing from loaded optimal IC for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    %solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                    IC = 'optimized';
                                    [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC,originalIC);
                                    disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                catch
                                    disp(['No saved solution. Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    tic
                                    [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                    disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                end
                            case 'off'
                                disp(['Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                tic
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                disp([num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's elapsed'])
                        end
                end
                pause(1)
        
                %%% save/inspect solution %%%
                switch run 
                    case 'energygrowth'
                        if exist('l2norms_avg','var') == 0
                            l2norms_avg = NaN(length(L1_scale),length(test_parameter)+1);
                            l2norms_mode = NaN(length(L1_scale),length(test_parameter)+1);
                            l2norms_guess = NaN(T/dt,numberoftests);
                            l2norms_opt = NaN(T/dt,numberoftests);
                        end
                        %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, K, L_s1, L_s2, 0);                   % load solution
                        %[u_energyL2, ~] = load_2DKSsolution('energyL2', IC, dt, T, N, K, L_s1, L_s2, 0, 0);                       % load solution
                        %l2norms_guess(1:size(u_energyL2,1),testcounter) = u_energyL2;
                        [u_energyL2, ~] = load_2DKSsolution('energyL2', 'optimized', dt, T, N, K, L_s1, L_s2, [tol T], IC);                       % load solution
                        l2norms_opt(1:size(u_energyL2,1),testcounter) = u_energyL2;
                        u_energyL2rd = round(u_energyL2,1);
                        l2norms_mode(domain_i,param_i) = mode(u_energyL2rd);
                        l2norms_avg(domain_i,param_i) = mean(u_energyL2(ceil(end/2):end,1));
                        if param_i == length(test_parameter) 
                            l2norms_avg(domain_i,param_i+1) = std(l2norms_avg(domain_i,1:param_i));
                            l2norms_mode(domain_i,param_i+1) = std(l2norms_mode(domain_i,1:param_i));
                        end
                    case {'L','N','dt','IC'}                
                        %save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, K, L_s1, L_s2,[0 T],0);                      % save solution
                        %plot_2DKS(save_each, 'state', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                                % save/inspect state                             % save/inspect terminal state
                        plot_2DKS(save_each, 'diagnostics', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                            % save/inspect dynamical characteristics
                        %plot_2DKS(save_each, 'norms', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                                  % save/inspect dynamical characteristics
                        plot_2DKS(save_each, 'gif', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,100);                                     % save/inspect time evolution 
                        disp([num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's elapsed'])
                        close all                                                                                                         % close any open figures
                    case 'kappa' 
                        %pertIC = IC;
                        if testcounter > 1
                            [~,~,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,kappalist,rieszlist,Ntime_save_max);
                        else
                            [~,~,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,0,0,Ntime_save_max);
                        end               
                        delete_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime_save_max T],0);
                        delete_2DKSsolution('backward', IC, dt, T, N, K, L_s1, L_s2, [Ntime_save_max T],0);
                    case 'optimize' 
                        [J_opt, J_history , v_TC_opt , u_IC_opt] = optimize_2DKS(optmethod,IC,N,K,L_s1,L_s2,dt,T,u_TC,v_TC,u_IC,Ntime_save_max,originalIC,tol);
                        switch continuation
                            case 'off'
                                if (J_opt - J_history(1))/J_history(1) < tol
                                    optmethod = 'RG';
                                    disp('RCG did not work, trying RG method...')
                                    [J_opt_RG, J_history_RG , v_TC_opt_RG , u_IC_opt_RG] = optimize_2DKS(optmethod,IC,N,K,L_s1,L_s2,dt,T,u_TC,v_TC,u_IC,Ntime_save_max,originalIC,tol);
                                    if J_opt_RG > J_opt
                                        J_opt = J_opt_RG;
                                        J_history = J_history_RG;
                                        v_TC_opt = v_TC_opt_RG;
                                        u_IC_opt = u_IC_opt_RG;
                                        RCGon = 0;
                                    end
                                    optmethod = 'RCG';
                                end
                        end
                        %u_pert = 'noise';
                        %[ ~ , ~ , LLE ] = solve_2DKS(IC,'LLE',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC_opt,u_pert);
                        % hessapprox = localoptimalitytest(J_opt,u_IC_opt,IC,N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max);
                        IC = strjoin(initialcondition(param_i),'');
                        %[ ~ , ~ , LLE_og ] = solve_2DKS(IC,'LLE',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC_opt,u_pert);
                        disp(['Solved optimization problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                        fprintf('Iterations \t Initial J \t Optimal J \t Wall Clock \n ')
                        fprintf('----------------------------------------------------------------------------------------------------------\n')
                        fprintf('%02d \t %.4f\t %.4f \t ', length(J_history), J_history(1,1), J_history(end,1) )
                        disp(datetime)
                        currentrow = testrow + (param_i-1)*length(initialKmagnitude)*length(timewindow)*length(L1_scale)*length(L2_scale); 
                        Jinitdata(currentrow,1) = K;
                        Joptdata(currentrow,1) = Jinitdata(currentrow,1);
                        Jinitdata(currentrow,2) = L_s1;
                        Joptdata(currentrow,2) = Jinitdata(currentrow,2);
                        Jinitdata(currentrow,3) = T;
                        Joptdata(currentrow,3) = Jinitdata(currentrow,3);
                        Jinitdata(currentrow,3) = J_history(1,1);
                        Joptdata(currentrow,4) = J_history(end,1);
                        save_2DKSsolution('optimal', u_IC_opt, v_TC_opt, 0, IC, dt, T, N, K, L_s1, L_s2, [1 T], tol); % save solution to machine
                        newopt = 1;
                        
                        maxL2inT = plot_2DKS(save_each, 'norms', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,u_IC_opt,u_IC,[tol,RCGon,100,newopt]);
                        Joptdata(currentrow,5:6) = maxL2inT;
                        save_measures('optimization', Jinitdata(currentrow,:), Joptdata(currentrow,:), 1, IC, N, dt, timewindow, K, L_s1, L_s2);
                        fprintf('Saved objective data to file at ')
                        fprintf('%01dh%02dm%02ds\t\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
                        if optfigs == 1
                            plot_2DKS(save_each, 'optdiag', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,u_IC_opt,u_IC,[tol,RCGon,100,newopt]);                                         
                        end
                        close all
                        switch optmethod
                            case 'RCG'
                                RCGon = 1;                              % RCG switch for media files
                            otherwise
                                RCGon = 0;
                        end
                        fprintf('Deleting all evolution data stored in directory...\n')
                        delete_2DKSsolution('forward', 'optimized', dt, T, N, K, L_s1, L_s2, [Ntime_save_max T],originalIC);
                        delete_2DKSsolution('backward', 'optimized', dt, T, N, K, L_s1, L_s2, [Ntime_save_max T],originalIC);
                        delete_2DKSsolution('backward', originalIC, dt, T, N, K, L_s1, L_s2, [Ntime_save_max T],originalIC);
                        fprintf('Optimization run complete.\n')
                end
            end
        end
        switch run 
            case 'kappa'    % designed for 5 or 10 tests only
                %plot_measures('kappa', dt, pertIC, N, timewindow, K, L_s1, L_s2, testcounter, length(timewindow));
            case 'optimize'
                %rowstart = testcounter/(length(initialcondition)) - length(timewindow) + 1;
                %rowend = testcounter/(length(initialcondition));
                %save_measures('optimization', Jinitdata(rowstart:rowend,:), Joptdata(rowstart:rowend,:), 1, initialcondition, N, dt, timewindow, K, L_s1, L_s2);
        end
        %% end window block
    end
    switch run 
        case 'optimize'
            if length(timewindow) > 1
                %plot_measures('optimization', dt, initialcondition, N, T, K, L_s1, L_s2, Jinitdata, Joptdata, IC, tol);
            end
            %save_measures('optimization', Jinitdata, Joptdata, numberoftests, initialcondition, N, dt, timewindow, K, L_scale, L_s2);
    end
end

switch run 
    case 'energygrowth'     % designed for comparing two initial conditions only 
        if exist('L_target','var')
            plot_measures('energygrowth', L1_scale, initialcondition, l2norms_mode, T, K, L_target, L_s2, l2norms_avg, length(timewindow), 0, 0);
            save_measures('energygrowth', l2norms_mode, l2norms_avg, 0, initialcondition, 0, K, L_scale, T, L_target, 0);
        end
    case 'N'                % spatial convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('spatial', dt, IC, gridsize, T, K, L_s1, L_s2, L_s2, length(timewindow), 0, 0);
        save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, K, L_s1, L_s2);
    case 'dt'               % temporal convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('temporal', timestep, IC, N, T, K, L_s1, L_s2, L_s2, length(timewindow), 0, 0);
        save_measures('temporal', error_2, error_inf, comptime, IC, N, 0, T, K, L_s1, L_s2);
end

end

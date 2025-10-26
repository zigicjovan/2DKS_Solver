<<<<<<< HEAD
%function main_2DKS
tic

%%% choose test settings %%%
run = 'IC';                               % switch to 'optimize', 'plotOptIC', 'energygrowth', 'L', 'N', 'dt', 'IC', 'kappa'
restart = 1;                                    % binary switch to generate new test counters
optfigs = 1;                                    % generate optimization diagnostic figures 
continuation = 'optIC';                           % 'IC' for optimal IC from file, 'forward' for optimized forward solution, 'off' to generate new data
optmethod = 'RCG';                              % RCG, RG, or RCGd5 (start after 5th iter)
Ntime_save_max = 10000;                         % choose maximum number of samples per file
timestep = 1e-4;                               % time-step sizes
gridsize = 64;                                  % grid sizes
tol = 1e-8;                                    % set optimization tolerance critera
optT = 10^(-1.5);

%%% choose parameter testing ranges %%%
initialKmagnitude = 10^(4);                        % initial L^2 energy magnitudes
L_scale = 1.1:.02:1.1;                        % domain sizes
timewindow = logspace(2.0,2.0,1);         % time windows
=======
function main_2DKS(dtc,Nc,Kstart,Kend,Knum,ellstart,ellend,ellgap,Tstart,Tend,Tnum)

tic

%%% choose test settings %%%
run = 'optimize';                               % switch to 'optimize', 'plotOptIC', 'energygrowth', 'L', 'N', 'dt', 'IC', 'kappa'
restart = 1;                                    % binary switch to generate new test counters
optfigs = 1;                                    % generate optimization diagnostic figures 
continuation = 'IC';                            % 'IC' for optimal IC from file, 'forward' for optimized forward solution, 'off' to generate new data
optmethod = 'RCG';                              % RCG, RG, or RCGd5 (start after 5th iter)
Ntime_save_max = 10000;                         % choose maximum number of samples per data file
timestep = dtc;                                 % time-step sizes
gridsize = Nc;                                  % grid sizes
tol = 1e-1;                                     % set optimization tolerance critera
optT = 10^(-1.5);                               % T parameter of optimal IC for asymptotic simulations

%%% choose parameter testing ranges %%%
initialKmagnitude = logspace(Kstart,Kend,Knum);                        % initial L^2 energy magnitudes
L_scale = ellstart:ellgap:ellend;                        % domain sizes
timewindow = logspace(Tstart,Tend,Tnum);         % time windows
>>>>>>> a5f484f (Initial upload)
initialcondition = {'s1'}; % initial conditions

%timewindow = timewindow(2);
switch run
    case 'kappa'
        kappapert = {'stg30'};                          % perturbation functions
        L_target = [2.36,2.78];                         % domain sizes of interest
end

%%% choose default parameters %%%
K = initialKmagnitude(1);                       % magnitude of initial condition L^2 energy, use K=0 for default IC norm
L_s1 = L_scale(1);                              % length-scale parameter in dim 1
L_s2 = L_scale(1);                              % length-scale parameter in dim 2
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
    numberoftests = length(initialcondition)*length(initialKmagnitude)*length(timewindow)*length(L_scale);
    testcounter = 0;
    testrow = 0;
    
    Joptdata = NaN(length(initialKmagnitude)*length(timewindow)*length(L_scale),3+3*length(initialcondition));
    Jinitdata = NaN(length(initialKmagnitude)*length(timewindow)*length(L_scale),3+1*length(initialcondition));
end

for energy_i = 1 : length(initialKmagnitude)

    switch run
        case 'kappa'
            pertIC = strjoin(kappapert(energy_i),'');
    end
    K = initialKmagnitude(energy_i);                 % adjust initial L^2 energy magnitude

    for domain_i = 1 : length(L_scale)
    
        L_s1 = L_scale(domain_i);                         % adjust length-scale parameter in dim 1
        L_s2 = L_s1;                                    % adjust length-scale parameter in dim 2
    
        %% start window block
        for window_i = 1 : length(timewindow)
    
            T = timewindow(window_i);                    % adjust simulation time window
    
            switch run                                  % set which test to run
                case 'L'
                    test_parameter = L_scale;           % (dynamical behavior)
                case 'N'
                    test_parameter = gridsize;          % (spatial convergence)
                case 'dt'
                    test_parameter = timestep;          % (temporal convergence)
                case {'IC', 'energygrowth','kappa','optimize','plotOptIC'}
                    test_parameter = initialcondition;  
            end
        
            testrow = testrow + 1;

            for param_i = 1 : length(test_parameter)          % length('x') indicates 'x' testing
        
                if window_i > length(timewindow)/2
                    Nfull = N;
                    dtfull = dt;
                    %N = 2/3*N;
                    %dt = 2*dt;
                end
                %%% set variable parameters %%%
                switch run 
                    case 'L'
                        L_s2 = L_scale(param_i);                      % length-scale parameter in dim 2 
                    case 'N'
                        N = gridsize(param_i);                        % number of grid points 
                    case 'dt'
                        dt = timestep(param_i);                       % length of time-step
                    case {'IC', 'energygrowth','optimize'}
                        IC = strjoin(initialcondition(param_i),'');   % choice of initial condition
                    case {'kappa'}
                        save_each = 1;
                        IC = strjoin(initialcondition(param_i),'');   % choice of initial condition
                end
                
                testcounter = testcounter + 1;
                disp(['Test ' num2str(testcounter) ' of ' num2str(numberoftests) ])
                %%% solve forward-time PDE problem %%%
                switch run 
                    case {'L','N','dt','IC','kappa','plotOptIC'} 
                        switch continuation
                            case 'optIC' % choose which parameters define optimized initial data
                                initIC = IC;
                                optdt = dt;
                                %optT = optT;
                                optN = N;
                                optK = K;
                                optL_s1 = L_s1;
                                optL_s2 = L_s2;
                                opttol = tol; 
                                [ u_IC , ~ ] = load_2DKSsolution('optimal', initIC, optdt, optT, optN, optK, optL_s1, optL_s2, opttol, 0);
                                tic
                                disp(['Using chosen optimal IC for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                IC = 'optimized';
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC,originalIC);
                                time = toc;
                                disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                            case 'forward'
                                try 
                                    time1 = ceil(T/(dt*save_each));
                                    if time1 > Ntime_save_max 
                                        timeIC = Ntime_save_max;
                                        timeTC = mod(time1,Ntime_save_max);
                                        time2 = ceil(Ntime_save_max*dt*save_each);
                                    else
                                        timeIC = time1;
                                        timeTC = time1;
                                        time2 = T;
                                    end
                                    [ u , ~ ] = load_2DKSsolution('forward', 'optimized', dt, time2, N, K, L_s1, L_s2, timeIC, IC);
                                    u_IC = u(:,1);
                                    tic
                                    disp(['Continuing from loaded optimized solution for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                catch
                                    disp(['No saved solution. Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    tic
                                    [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                    time = toc;
                                    disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                end
                            case 'off'
                                disp(['Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                tic
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                time = toc;
                                disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                        end
                    case 'optimize'
                        originalIC = IC;
                        switch continuation
                            case 'IC'
                                try
                                    [ u_IC , ~ ] = load_2DKSsolution('optimal', IC, dt, T, N, K, L_s1, L_s2, tol, 0);
                                    tic
                                    disp(['Continuing from loaded optimal IC for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                    IC = 'optimized';
                                    [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC,originalIC);
                                    time = toc;
                                    disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                catch
                                    try
                                        time1 = ceil(T/(dt*save_each));
                                        if time1 > Ntime_save_max 
                                            timeIC = Ntime_save_max;
                                            timeTC = mod(time1,Ntime_save_max);
                                            time2 = ceil(Ntime_save_max*dt*save_each);
                                        else
                                            timeIC = time1;
                                            timeTC = time1;
                                            time2 = T;
                                        end
                                        [ u , ~ ] = load_2DKSsolution('forward', 'optimized', dt, time2, N, K, L_s1, L_s2, timeIC, IC);
                                        [ ~ , v ] = load_2DKSsolution('forward', 'optimized', dt, T, N, K, L_s1, L_s2, timeTC, IC);
                                        u_IC = u(:,1);
                                        v_TC = v(:,end);
                                        u_TC = real(ifft2(v_TC));
                                        IC = 'optimized';
                                        tic
                                        disp(['No optimal IC file found. Continuing from loaded optimized solution for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    catch
                                        disp(['No saved solution. Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                        tic
                                        [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                        time = toc;
                                        disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                    end
                                end
                            case 'forward'
                                try
                                    time1 = ceil(T/(dt*save_each));
                                    if time1 > Ntime_save_max 
                                        timeIC = Ntime_save_max;
                                        timeTC = mod(time1,Ntime_save_max);
                                        time2 = ceil(Ntime_save_max*dt*save_each);
                                    else
                                        timeIC = time1;
                                        timeTC = time1;
                                        time2 = T;
                                    end
                                    [ u , ~ ] = load_2DKSsolution('forward', 'optimized', dt, time2, N, K, L_s1, L_s2, timeIC, IC);
                                    [ ~ , v ] = load_2DKSsolution('forward', 'optimized', dt, T, N, K, L_s1, L_s2, timeTC, IC);
                                    u_IC = u(:,1);
                                    v_TC = v(:,end);
                                    u_TC = real(ifft2(v_TC));
                                    IC = 'optimized';
                                    tic
                                    disp(['Continuing from loaded optimized solution for forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                catch
                                    disp(['No saved solution. Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                    tic
                                    [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                    time = toc;
                                    disp(['Solved forward problem at ' num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's'])
                                end
                            case 'off'
                                disp(['Solving forward-time problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                                tic
                                [ v_TC , u_TC , u_IC ] = solve_2DKS(IC,'forward',N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,0,0);
                                time = toc;
                                disp([num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's elapsed'])
                        end
                end
                pause(1)
        
                %%% save/inspect solution %%%
                switch run 
                    case 'energygrowth'
                        if exist('l2norms_avg','var') == 0
                            l2norms_avg = NaN(length(L_scale),length(test_parameter)+1);
                            l2norms_mode = NaN(length(L_scale),length(test_parameter)+1);
                            l2norms_guess = NaN(T/dt,numberoftests);
                            l2norms_opt = NaN(T/dt,numberoftests);
                        end
                        %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, K, L_s1, L_s2, 0);                   % load solution
                        %[u_energyL2, ~] = load_2DKSsolution('energyL2', IC, dt, T, N, K, L_s1, L_s2, 0, 0);                       % load solution
                        %l2norms_guess(1:size(u_energyL2,1),testcounter) = u_energyL2;
                        [u_energyL2, ~] = load_2DKSsolution('energyL2', 'optimized', dt, T, N, K, L_s1, L_s2, tol, IC);                       % load solution
                        l2norms_opt(1:size(u_energyL2,1),testcounter) = u_energyL2;
                        u_energyL2rd = round(u_energyL2,1);
                        l2norms_mode(domain_i,param_i) = mode(u_energyL2rd);
                        l2norms_avg(domain_i,param_i) = mean(u_energyL2(ceil(end/2):end,1));
                        if param_i == length(test_parameter) 
                            l2norms_avg(domain_i,param_i+1) = std(l2norms_avg(domain_i,1:param_i));
                            l2norms_mode(domain_i,param_i+1) = std(l2norms_mode(domain_i,1:param_i));
                        end
                    case {'L','N','dt','IC'}                
                        %save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, K, L_s1, L_s2,0,0);                      % save solution
                        %plot_2DKS(save_each, 'initial', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                                % save/inspect initial state
                        %plot_2DKS(save_each, 'terminal', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                               % save/inspect terminal state
                        plot_2DKS(save_each, 'diagnostics', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                            % save/inspect dynamical characteristics
                        %plot_2DKS(save_each, 'norms', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,0);                                  % save/inspect dynamical characteristics
                        plot_2DKS(save_each, 'gif', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max, 0,100);                                     % save/inspect time evolution 
                        disp([num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's elapsed'])
                        close all                                                                                                         % close any open figures
                    case {'plotOptIC'}                
                        plot_2DKS(save_each, 'optdiag', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,[tol,RCGon,100]); 
                        %{
                        plot_2DKS(save_each, 'diagnostics', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,tol); 
                        plot_2DKS(save_each, 'gif', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,100);
                        plot_2DKS(save_each, 'gif', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,100); 
                        plot_2DKS(save_each, 'initial', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,0);
                        plot_2DKS(save_each, 'terminal', IC, N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,0);
                        plot_2DKS(save_each, 'initial', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,tol);
                        plot_2DKS(save_each, 'terminal', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,tol);  
                        %}
                        close all
                    case 'kappa' 
                        %pertIC = IC;
                        if testcounter > 1
<<<<<<< HEAD
                            [kappa,gat_riesz,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,kappalist,rieszlist,Ntime_save_max);
                        else
                            [kappa,gat_riesz,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,0,0,Ntime_save_max);
=======
                            [~,~,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,kappalist,rieszlist,Ntime_save_max);
                        else
                            [~,~,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,0,0,Ntime_save_max);
>>>>>>> a5f484f (Initial upload)
                        end               
                        delete_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_save_max,0);
                        delete_2DKSsolution('backward', IC, dt, T, N, K, L_s1, L_s2, Ntime_save_max,0);
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
                        IC = strjoin(initialcondition(param_i),'');
<<<<<<< HEAD
                        disp([num2str(floor(toc/3600)) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's elapsed'])
                        disp(['Solved optimization problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                        disp(['Initial objective functional value: ' num2str(J_history(1,1))])
                        disp(['Optimal objective functional value: ' num2str(J_history(end,1))])
                        disp(['Number of iterations: ' num2str(length(J_history))])
=======
                        disp(['Solved optimization problem for K = ' num2str(K) ', L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt) ', N = ' num2str(N)])
                        fprintf('Iterations \t Initial J \t Optimal J \t Wall Clock \n ')
                        fprintf('----------------------------------------------------------------------------------------------------------\n')
                        fprintf('%02d \t %.4f\t %.4f \t ', length(J_history), J_history(1,1), J_history(end,1) )
>>>>>>> a5f484f (Initial upload)
                        disp(datetime)
                        Jinitdata(testrow,1) = K;
                        Joptdata(testrow,1) = Jinitdata(testrow,1);
                        Jinitdata(testrow,2) = L_s1;
                        Joptdata(testrow,2) = Jinitdata(testrow,2);
                        Jinitdata(testrow,3) = T;
                        Joptdata(testrow,3) = Jinitdata(testrow,3);
                        Jinitdata(testrow,param_i+3) = J_history(1,1);
                        Joptdata(testrow,(param_i-1)*3+4) = J_history(end,1);
                        save_2DKSsolution('optimal', u_IC_opt, v_TC_opt, 0, IC, dt, T, N, K, L_s1, L_s2, 1, tol); % save solution to machine
<<<<<<< HEAD
=======
                        fprintf('Saving diagnostic figures... ')
>>>>>>> a5f484f (Initial upload)
                        if optfigs == 1
                            maxL2inT = plot_2DKS(save_each, 'optdiag', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,[tol,RCGon,100]);                       
                        else
                            maxL2inT = plot_2DKS(save_each, 'norms', 'optimized', N, dt, T, K, L_s1, L_s2,Ntime_save_max,IC,[tol,RCGon,100]);                     
                        end
                        Joptdata(testrow,(param_i-1)*3+5:(param_i-1)*3+6) = maxL2inT;  
                        close all
                        switch optmethod
                            case 'RCG'
                                RCGon = 1;                              % RCG switch for media files
                            otherwise
                                RCGon = 0;
                        end
                        if window_i > length(timewindow)/2
                            N = Nfull;
                            dt = dtfull;
                        end
<<<<<<< HEAD
=======
                        fprintf('Optimization run complete.')
>>>>>>> a5f484f (Initial upload)
                end
            end
        end
        switch run 
            case 'kappa'    % designed for 5 or 10 tests only
                plot_measures('kappa', dt, pertIC, N, timewindow, K, L_s1, L_s2, testcounter, length(timewindow));
            case 'optimize'
                rowstart = testcounter/(length(initialcondition)) - length(timewindow) + 1;
                rowend = testcounter/(length(initialcondition));
                save_measures('optimization', Jinitdata(rowstart:rowend,:), Joptdata(rowstart:rowend,:), 1, initialcondition, N, dt, timewindow, K, L_s1, L_s2);
        end
        %% end window block
    end
    switch run 
        case 'optimize'
            if length(timewindow) > 1
                plot_measures('optimization', dt, initialcondition, N, T, K, L_s1, L_s2, Jinitdata, Joptdata, IC, tol);
            end
            save_measures('optimization', Jinitdata, Joptdata, numberoftests, initialcondition, N, dt, timewindow, K, L_scale, L_s2);
    end
end

switch run 
    case 'energygrowth'     % designed for comparing two initial conditions only 
        if exist('L_target','var')
            plot_measures('energygrowth', L_scale, initialcondition, l2norms_mode, T, K, L_target, L_s2, l2norms_avg, length(timewindow), 0, 0);
            save_measures('energygrowth', l2norms_mode, l2norms_avg, 0, initialcondition, 0, K, L_scale, T, L_target, 0);
        end
    case 'N'                % spatial convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('spatial', dt, IC, gridsize, T, K, L_s1, L_s2, L_s2, length(timewindow), 0, 0);
        save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, K, L_s1, L_s2);
    case 'dt'               % temporal convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('temporal', timestep, IC, N, T, K, L_s1, L_s2, L_s2, length(timewindow), 0, 0);
        save_measures('temporal', error_2, error_inf, comptime, IC, N, 0, T, K, L_s1, L_s2);
<<<<<<< HEAD
=======
end

>>>>>>> a5f484f (Initial upload)
end

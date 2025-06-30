%function main_2DKS
tic

%%% choose test %%%
run = 'kappa';                                     % switch to 'L', 'N', 'dt', 'IC', 'kappa', 'energygrowth'

%%% choose parameter testing ranges %%%
L_scale =  2.36 ;        % domain sizes
timestep = [ .005, .001, .0005 , .0001 ];                                % time-step sizes
gridsize = 48;                                  % grid sizes
timewindow = linspace(50,500,10);          % time windows
initialcondition = { 's1' };                  % initial conditions
kappapert = 0;                  % perturbation function
L_target = 2.36;                   % domain size of interest

%%% choose default parameters %%%
L_s1 = L_scale(1);                              % length-scale parameter in dim 1
L_s2 = L_scale(1);                              % length-scale parameter in dim 2
dt = timestep(4);                               % length of time-step
N = gridsize(1);                                % number of grid points
T = timewindow(1);                              % length of simulation time window
IC = strjoin(initialcondition(1),'');           % initial condition
save_each = 1/dt;                               % number of iterations between saved timepoints - 1/dt to save each 1 T

numberoftests = length(initialcondition)*length(kappapert)*length(timewindow)*length(L_scale);
testcounter = 0;

%for init = 4 : length(kappapert)

    %pertIC = strjoin(kappapert(init),'');

    for choros = 1 : length(L_scale)
    
        L_s1 = L_scale(choros);                     % adjust length-scale parameter in dim 1
        L_s2 = L_s1;                                % adjust length-scale parameter in dim 2
    
        for chronos = 1 : length(timewindow)
    
            T = timewindow(chronos);                % adjust simulation time window
    
            switch run                                  % set which test to run
                case 'L'
                    test_parameter = L_scale;           % (dynamical behavior)
                case 'N'
                    test_parameter = gridsize;          % (spatial convergence)
                case 'dt'
                    test_parameter = timestep;          % (temporal convergence)
                case {'IC', 'energygrowth'}
                    test_parameter = initialcondition;  
                case 'kappa'
                    test_parameter = initialcondition; 
                    epscheck = logspace(-15,-1,15);
                    J_pert = NaN(length(epscheck),1);                 % store perturbed objective functionals  
                    gateaux_deriv = NaN(length(epscheck),1);          % store kappa test numerators 
                    kappa = NaN(length(epscheck),1);                  % store kappa test results 
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
                    case 'kappa'
                        save_each = 1;                          % save all timesteps for backward solver
                        IC = strjoin(initialcondition(k),'');   % choice of initial condition
                        pertIC = IC;
                end
                
                testcounter = testcounter + 1;
                disp(['Test ' num2str(testcounter) ' of ' num2str(numberoftests) ])
                %%% solve forward-time PDE problem %%%
                switch run 
                    case {'L','N','dt','IC','kappa'} 
                        disp(['Solving forward-time problem for L = ' num2str(L_s1) ', T = ' num2str(T) ', IC = ' IC ', dt = ' num2str(dt)])
                        tic
                        [ v_TC , u_TC ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,0,0);
                        time = toc;
                        toc
                end
        
                %%% save/inspect solution %%%
                switch run 
                    case 'energygrowth'
                        if exist('l2norms_avg','var') == 0
                            l2norms_avg = NaN(length(L_scale),length(test_parameter)+1);
                            l2norms_mode = NaN(length(L_scale),length(test_parameter)+1);
                        end
                        %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2, 0);               % load solution
                        [u_normL2, ~] = load_2DKSsolution('normL2', IC, dt, T, N, L_s1, L_s2, 0);                   % load solution
                        u_normL2rd = round(u_normL2,1);
                        l2norms_mode(choros,k) = mode(u_normL2rd);
                        l2norms_avg(choros,k) = mean(u_normL2(ceil(end/2):end,1));
                        if k == length(test_parameter) 
                            l2norms_avg(choros,k+1) = std(l2norms_avg(choros,1:k));
                            l2norms_mode(choros,k+1) = std(l2norms_mode(choros,1:k));
                        end
                    case {'L','N','dt','IC'}                
                        %save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, L_s1, L_s2);               % save solution
                        %plot_2DKS(u_n, 'initial', IC, N, dt, T, L_s1, L_s2, 0,0);                               % save/inspect initial state
                        %plot_2DKS(u_n, 'terminal', IC, N, dt, T, L_s1, L_s2, 0,0);                              % save/inspect terminal state
                        %plot_2DKS(u_n, 'diagnostics', IC, N, dt, T, L_s1, L_s2, 0,0);                           % save/inspect dynamical characteristics
                        %plot_2DKS(u_n, 'norms', IC, N, dt, T, L_s1, L_s2, 0,0);                                 % save/inspect dynamical characteristics
                        plot_2DKS(save_each, 'gif', IC, N, dt, T, L_s1, L_s2, 0,0);                                    % save/inspect time evolution 
                        toc
                        close all                                                                                % close any open figures
                    case 'kappa' 
                        %plot_2DKS(save_each, 'gif', IC, N, dt, T, L_s1, L_s2, 0,0); 
                        if exist('kappalist','var') == 0
                            kappalist = NaN(length(epscheck),numberoftests); 
                            rieszlist = NaN(numberoftests,1); 
                        end
                        L1 = 2*pi*L_s1;
                        L2 = 2*pi*L_s2;
                        J_init = sum( u_TC .* conj(u_TC) )*(L1*L2)/N^2;                                          % initial objective functional (L^2 inner product of terminal forward state)        
                        disp(['Solving adjoint problem for pertIC = ' pertIC])
                        [v_adjIC, u_adjIC, u_pertIC] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,v_TC,pertIC);    % solve adjoint equation
                        toc
                        gat_riesz = sum( u_adjIC .* conj(u_pertIC) )*(L1*L2)/N^2;                               % kappa test denominator 
                        rieszlist(testcounter,1) = abs(gat_riesz);                                              % save kappa test denominator values
                        gradJ = v_adjIC;                                                                        % objective gradient
                        save_each = 1/dt;                                                                       % release memory
                        disp('Solving perturbed forward-time problems...')
                        for i = 1:length(epscheck)   
                            eps = epscheck(i);                                                                  % define perturbation 
                            [ ~ , u_pertTC ] = solve_2DKS(IC,'kappa',N,L_s1,L_s2,dt,T,save_each,eps,pertIC);    % solve perturbed forward equation
                            J_pert(i,1) = sum( u_pertTC .* conj(u_pertTC) )*(L1*L2)/N^2;                        % perturbed objective functional
                            gateaux_deriv(i,1) = (J_pert(i,1) - J_init)/eps;                                    % kappa test numerator 
                            kappa(i,1) = gateaux_deriv(i,1)/gat_riesz;                                          % kappa test numerator 
                            toc
                        end
                        plot_2DKS(save_each, 'kappa', IC, N, dt, T, L_s1, L_s2, kappa,pertIC);                        % save/inspect kappa test figure
                        kappalist(:,testcounter) = kappa;                                                       % save kappa test values
                        close all                                                                               % close any open figures
                        kappalist_file = [pwd '/data/kappa/kappalist_' IC '_p' pertIC '_N_' num2str(N) '' ...
                        '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                        writematrix(kappalist, kappalist_file,'Delimiter','tab');
                        rieszlist_file = [pwd '/data/kappa/rieszlist_' IC '_p' pertIC '_N_' num2str(N) '' ...
                         '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                        writematrix(rieszlist, rieszlist_file,'Delimiter','tab');
                end
            end
        end
    end
%end

switch run 
    case 'energygrowth'                                                                             % currently for comparing two initial conditions only 
        plot_measures('energygrowth', L_scale, initialcondition, l2norms_mode, T, L_target, l2norms_avg);
        save_measures('energygrowth', l2norms_mode, l2norms_avg, 0, initialcondition, 0, L_scale, T, L_target, 0);
    case 'N'                                                                                        % spatial convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('spatial', dt, IC, gridsize, T, L_s1, L_s2);
        save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, L_s1, L_s2);
    case 'dt'                                                                                       % temporal convergence: error analysis and computational time
        [error_2,error_inf,comptime] = plot_measures('temporal', timestep, IC, N, T, L_s1, L_s2);
        save_measures('temporal', error_2, error_inf, comptime, IC, N, 0, T, L_s1, L_s2);
end
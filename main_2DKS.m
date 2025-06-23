%function main_2DKS
tic

%%% choose test %%%
run = 'kappa';                                     % switch to 'L', 'N', 'dt', 'IC', 'kappa', 'load

%%% choose parameter testing ranges %%%
L_scale =  [ 2.36,2.78] ;        % domain sizes
timestep = [ .005, .001, .0005 ];                                % time-step sizes
gridsize = 48;                                  % grid sizes
timewindow = linspace(20,200,10);          % time windows
initialcondition = { 's1', 'stg1' };                  % initial conditions
kappapert = 0;%{ 's1', 'stg1', 'stg30' , 's30' };                  % perturbation function

%%% choose default parameters %%%
L_s1 = L_scale(1);                              % length-scale parameter in dim 1
L_s2 = L_scale(1);                              % length-scale parameter in dim 2
dt = timestep(2);                               % length of time-step
N = gridsize(1);                                % number of grid points
T = timewindow(1);                              % length of simulation time window
IC = strjoin(initialcondition(1),'');           % initial condition
save_each = 1/dt;                               % number of iterations between saved timepoints - 1/dt to save each 1 T

numberoftests = length(initialcondition)*length(kappapert)*length(timewindow)*length(L_scale);
kappalist = NaN(15,numberoftests); % temporary for kappa tests
rieszlist = NaN(numberoftests,1); % temporary for kappa tests
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
                case {'IC', 'load'}
                    test_parameter = initialcondition;  
                case 'kappa'
                    test_parameter = initialcondition; 
                    J_pert = NaN(15,1);                 % store perturbed objective functionals  
                    gateaux_deriv = NaN(15,1);          % store kappa test numerators 
                    kappa = NaN(15,1);                  % store kappa test results 
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
                    case {'IC', 'load'}
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
                        [ v_n , u_n ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,0,0);
                        time = toc;
                        toc
                end
        
                %%% save/inspect solution %%%
                switch run 
                    case 'load'
                        if exist('l2norms_avg','var') == 0
                            l2norms_avg = NaN(length(L_scale),length(test_parameter));
                            l2norms_max = NaN(length(L_scale),length(test_parameter));
                            l2norms_mode = NaN(length(L_scale),length(test_parameter));
                            l2norms_diff = NaN(length(L_scale),3);
                        end
                        %[u_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2);               % load solution
                        [u_normL2, ~] = load_2DKSsolution('normL2', IC, dt, T, N, L_s1, L_s2);                       % load solution
                        u_normL2rd = round(u_normL2,1);
                        l2norms_mode(choros,k) = mode(u_normL2rd);
                        l2norms_avg(choros,k) = mean(u_normL2(ceil(end/2):end,1));
                        l2norms_max(choros,k) = max(u_normL2);
                        if k == length(test_parameter) 
                            l2norms_diff(choros,1) = abs(l2norms_avg(choros,1) - l2norms_avg(choros,2));
                            l2norms_diff(choros,2) = abs(l2norms_max(choros,1) - l2norms_max(choros,2));
                            l2norms_diff(choros,3) = abs(l2norms_mode(choros,1) - l2norms_mode(choros,2));
                        end
                    case {'L','N','dt','IC'}                
                        %save_2DKSsolution('time_evolution', u_n, time, IC, dt, T, N, L_s1, L_s2);               % save solution
                        %plot_2DKS(u_n, 'initial', IC, N, dt, T, L_s1, L_s2, 0,0);                                 % save/inspect initial state
                        %plot_2DKS(u_n, 'terminal', IC, N, dt, T, L_s1, L_s2, 0,0);                                % save/inspect terminal state
                        %plot_2DKS(u_n, 'diagnostics', IC, N, dt, T, L_s1, L_s2, 0,0);                             % save/inspect dynamical characteristics
                        %plot_2DKS(u_n, 'norms', IC, N, dt, T, L_s1, L_s2, 0,0);                                   % save/inspect dynamical characteristics
                        plot_2DKS(u_n, 'gif', IC, N, dt, T, L_s1, L_s2, 0,0);                                     % save/inspect time evolution 
                        toc
                        close all                                                                                % close any open figures
                    case 'kappa' 
                        L1 = 2*pi*L_s1;
                        L2 = 2*pi*L_s2;
                        J_init = sum( u_n(:,end) .* conj(u_n(:,end)) )*(L1*L2)/N^2;                             % initial objective functional (L^2 inner product of terminal forward state)        
                        disp(['Solving adjoint problem for pertIC = ' pertIC])
                        [v_adj, u_adj, u_pertIC] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,v_n,pertIC);    % solve adjoint equation
                        toc
                        gat_riesz = sum( u_adj(:,1) .* conj(u_pertIC) )*(L1*L2)/N^2;                            % kappa test denominator 
                        rieszlist(testcounter,1) = gat_riesz;                               % temporary for kappa tests
                        gradJ = v_adj(:,1);                                                                     % objective gradient
                        save_each = 1/dt;                                                                       % release memory
                        clear v_n v_adj u_adj                                                                   % release memory
                        disp('Solving perturbed forward-time problems...')
                        for i = -15:1:-1
                            eps = 10^(i);                                                                       % define perturbation 
                            [ ~ , u_pert ] = solve_2DKS(IC,'kappa',N,L_s1,L_s2,dt,T,save_each,eps,pertIC);             % solve perturbed forward equation
                            J_pert(i+16,1) = sum( u_pert(:,end) .* conj(u_pert(:,end)) )*(L1*L2)/N^2;           % perturbed objective functional
                            gateaux_deriv(i+16,1) = (J_pert(i+16,1) - J_init)/eps;                              % kappa test numerator 
                            kappa(i+16,1) = gateaux_deriv(i+16,1)/gat_riesz;                                    % kappa test numerator 
                            toc
                        end
                        clear v_pert u_pert                                                                     % release memory
                        plot_2DKS(u_n, 'kappa', IC, N, dt, T, L_s1, L_s2, kappa,pertIC);                               % save/inspect kappa test figure
                        kappalist(:,testcounter) = kappa;                 % temporary for kappa tests
                        close all                                                                                % close any open figures
                end
            end
        end
    end
%end

switch run 
    case 'kappa'
        kappalist_file = [pwd '/data/kappa/kappalist_' IC '_p' pertIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        writematrix(kappalist, kappalist_file,'Delimiter','tab');
        rieszlist_file = [pwd '/data/kappa/rieszlist_' IC '_p' pertIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        writematrix(rieszlist, rieszlist_file,'Delimiter','tab');
    case 'load'
        %{
        figure
        plot(L_scale,l2norms_T(:,1))
        hold on
        plot(L_scale,l2norms_T(:,2))
        plot(L_scale,l2norms_T(:,3))
        plot(L_scale,l2norms_T(:,4))
        plot(L_scale,l2norms_T(:,5))
        plot(L_scale,l2norms_T(:,6))
        title("terminal")
        legend('s1', 's30', 'stg1', 'stg30' , 'tg1', 'tg30','Location','northwest')
        hold off

        figure
        plot(L_scale,l2norms_std(:,1))
        hold on
        plot(L_scale,l2norms_std(:,2))
        plot(L_scale,l2norms_std(:,3))
        title("std")
        legend('terminal', 'max', 'combined','Location','northwest')
        hold off
        %}

        figure
        plot(L_scale,l2norms_mode(:,1),'r-x')
        hold on
        plot(L_scale,l2norms_mode(:,2),'b-o')
        plot(L_scale,l2norms_avg(:,1),'g-+')
        plot(L_scale,l2norms_avg(:,2),'k-s')
        xl1 = xline(2.36,'--',{'2.36'});
        xl2 = xline(2.78,'--',{'2.78'});
        xl1.LabelVerticalAlignment = 'top';
        xl1.LabelHorizontalAlignment = 'left';
        xl2.LabelVerticalAlignment = 'top';
        xl2.LabelHorizontalAlignment = 'left';
        title('Modal energy ($t\in [0,500]$) and mean energy ($t\in [250,500]$)','Interpreter','latex')
        xlabel('Domain factor $\ell$','Interpreter','latex')
        %ylabel('$Energy value || \varphi_{IC}(t) ||_{L^2}$','Interpreter','latex')
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white') 
        legend('Modal $|| \varphi_{s1}(t) ||_{L^2}$',...
            'Modal $|| \varphi_{stg1}(t) ||_{L^2}$' ,...
            'Mean $|| \varphi_{s1}(t) ||_{L^2}$',...
            'Mean $|| \varphi_{stg1}(t) ||_{L^2}$' ,...
            'Location','northwest','Interpreter','latex')
        legend('boxoff')
        hold off

        figure
        plot(L_scale,l2norms_diff(:,1),'r-*')
        hold on
        plot(L_scale,l2norms_diff(:,3),'b-o')
        xl1 = xline(2.36,'--',{'2.36'});
        xl2 = xline(2.78,'--',{'2.78'});
        xl1.LabelVerticalAlignment = 'top';
        xl1.LabelHorizontalAlignment = 'left';
        xl2.LabelVerticalAlignment = 'top';
        xl2.LabelHorizontalAlignment = 'left';
        title('Difference in $|| \varphi_{IC}(t) ||_{L^2}$ values for ${IC} = \{s1,stg1 \}$ within $T =500$','Interpreter','latex')
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')  
        xlabel('Domain factor $\ell$','Interpreter','latex')
        legend('Difference in mean values','Difference in modal values','Location','northwest')

    case 'N'                                                                                        % spatial convergence: error analysis and computational time
        [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, IC, dt, T, L_s1, L_s2);
        save_measures('spatial', error_2, error_inf, comptime, IC, 0, dt, T, L_s1, L_s2);
    case 'dt'                                                                                       % temporal convergence: error analysis and computational time
        [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, IC, N, T, L_s1, L_s2);
        save_measures('temporal', error_2, error_inf, comptime, meICthod, N, 0, T, L_s1, L_s2);
end
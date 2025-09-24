function [kappa,gat_riesz,kappalist,rieszlist] = test_kappa(numberoftests,testcounter,IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,pertIC,kappalist,rieszlist,Ntime_save_max)

    if not(isfolder([pwd  '/data/kappa' ]))                                           % create local directories for data storage
        mkdir([pwd  '/data/kappa' ]);
        addpath([pwd  '/data/kappa' ]);
    end

    save_each = 1;                                                                          % save all timesteps for backward solver
    epscheck = logspace(-15,-1,15);
    J_pert = NaN(length(epscheck),1);                                                       % store perturbed objective functionals  
    gateaux_deriv = NaN(length(epscheck),1);                                                % store kappa test numerators 
    kappa = NaN(length(epscheck),1);                                                        % store kappa test results 
    if isscalar(kappalist)
        kappalist = NaN(length(epscheck),numberoftests);                                    % store multiple kappa tests
        rieszlist = NaN(numberoftests,1);                                                   % store kappa denominator values
    end
    L1 = 2*pi*L_s1;
    L2 = 2*pi*L_s2;
    J_init = sum( u_TC .* conj(u_TC) )*(L1*L2)/N^2;                                          % initial objective functional (L^2 inner product of terminal forward state)        
    disp(['Solving adjoint problem for pertIC = ' pertIC])
    [~, u_adjIC, u_pertIC] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,Ntime_save_max,v_TC,pertIC);    % solve adjoint equation
    toc
    gat_riesz = sum( u_adjIC .* conj(u_pertIC) )*(L1*L2)/N^2;                               % kappa test denominator 
    rieszlist(testcounter,1) = abs(gat_riesz);                                              % save kappa test denominator values
    save_each = 1/dt;                                                                       % release memory if necessary
    disp('Solving perturbed forward-time problems...')
    for i = 1:length(epscheck)   
        eps = epscheck(i);                                                                  % define perturbation 
        [ ~ , u_pertTC ] = solve_2DKS(IC,'kappa',N,L_s1,L_s2,dt,T,save_each,Ntime_save_max,eps,pertIC);    % solve perturbed forward equation
        J_pert(i,1) = sum( u_pertTC .* conj(u_pertTC) )*(L1*L2)/N^2;                        % perturbed objective functional
        gateaux_deriv(i,1) = (J_pert(i,1) - J_init)/eps;                                    % kappa test numerator 
        kappa(i,1) = gateaux_deriv(i,1)/gat_riesz;                                          % kappa test numerator 
        toc
    end
    plot_2DKS(save_each, 'kappa', IC, N, dt, T, L_s1, L_s2,Ntime_save_max,kappa,pertIC);    % save/inspect kappa test figure
    kappalist(:,testcounter) = kappa;                                                       % save kappa test values
    close all                                                                               % close any open figures
    kappalist_file = [pwd '/data/kappa/kappalist_' IC '_p' pertIC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(kappalist, kappalist_file,'Delimiter','tab');
    rieszlist_file = [pwd '/data/kappa/rieszlist_' IC '_p' pertIC '_N_' num2str(N) '' ...
     '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(rieszlist, rieszlist_file,'Delimiter','tab');

end
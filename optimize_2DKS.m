function [J_cur , J_history , u_TC , u_IC] = optimize_2DKS(IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,u_IC)

    maxiter = 1000;                                 % set maximum number of iterations
    if exist('kappalist','var') == 0
        %epscheck = logspace(-15,-1,15);
        J_history = NaN(maxiter,1);                        % store updated objective functionals  
        J_change = NaN(maxiter,1);                         % store changes in objective functional value
        stepsize_history = NaN(maxiter,1);                 % store optimization step sizes
        manifold_history = NaN(maxiter,1);                 % store manifold sizes
        %gateaux_deriv = NaN(length(epscheck),1);          % store kappa test numerators 
        %kappa = NaN(length(epscheck),1);                  % store kappa test results 
        %Jlist = NaN(length(epscheck),numberoftests); 
        %rieszlist = NaN(numberoftests,1); 
    end

    save_each = 1;                                                                          % save all timesteps for backward solver
    tol = 1e-5;                                                                             % set optimization tolerance critera
    iter = 1;                                                                               % start counting iterations  
    L1 = 2*pi*L_s1;                                                                         % dimension 1 length
    L2 = 2*pi*L_s2;                                                                         % dimension 2 length

    manifold_size = sum( u_IC .* conj(u_IC) )*(L1*L2)/N^2;                                  % current manifold (L^2 inner product of initial forward state)        
    J_cur = sum( u_TC .* conj(u_TC) )*(L1*L2)/N^2;                                          % current objective functional (L^2 inner product of terminal forward state) 
    J_history(iter,1) = J_cur;                                                              % store objective functional history
    manifold_history(iter,1) = manifold_size;                                               % store manifold sizes
         
    dir_old = 0;                                                                            % initialize old direction
    vectransport = 0;                                                                       % initialize current vector transport operator
    momentum_size = 0;                                                                      % initialize current momentum
    projGradJ_old = 0;                                                                      % initialize old projected objective gradient
    J_change(iter,1) = 1;                                                                   % initialize change in objective functional value 
    
    while (abs(J_change(iter,1)) > tol) && (iter <= maxiter)

        disp(['Solving adjoint problem for iteration = ' num2str(iter)])
        [~, GradJ_cur] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,v_TC,0);                   % current objective gradient via adjoint equation
        toc
        angleGradJ_cur = sum( u_IC .* conj(GradJ_cur) )*(L1*L2)/N^2;                                    % angle with current objective gradient
        proj_cur = GradJ_cur - (angleGradJ_cur/manifold_size).*(u_IC);                                  % current projection operator
        projGradJ_cur = proj_cur .* GradJ_cur;                                                          % current projected objective gradient
        if iter > 1
            angleDir_old = sum( u_IC .* conj(angleDir_old) )*(L1*L2)/N^2;                               % angle with old direction
            vectransport = (dir_old - (angleDir_old/manifold_size).*(u_IC))/manifold_size;              % current vector transport operator
            diff_projGradJ = projGradJ_cur - ( vectransport .* projGradJ_old );                         % momentum parameter term
            diff_projGradJ_ip = sum( projGradJ_cur .* conj(diff_projGradJ) )*(L1*L2)/N^2;               % momentum parameter numerator
            projGradJ_old_ip = sum( projGradJ_old .* conj(projGradJ_old) )*(L1*L2)/N^2;                 % momentum parameter denominator
            momentum_size = diff_projGradJ_ip / projGradJ_old_ip;                                       % current momentum parameter
        else
            GradJ_cur_ip = sum( GradJ_cur .* conj(GradJ_cur) )*(L1*L2)/N^2;                             % inner product current objective gradient
            step_size = .1;%-2*real(angleGradJ_cur/GradJ_cur_ip);                                        % initialize current step-size
            stepsize_history(iter,1) = step_size;                                                       % store optimization step size
        end
        dir_cur = projGradJ_cur + (momentum_size .* vectransport .* dir_old);                           % current direction
        [step_size,iter_search,J_search] = optimize_stepsize(dir_cur,u_IC,step_size,IC,N,L_s1,L_s2,dt,T); % current step-size via Brent's method
        disp(['Line-search problem number of iterations = ' num2str(iter_search)])
        J_history(iter+1:iter + iter_search,1) = J_search;                                              % store objective functional history from line search
        iter = iter + iter_search;                                                                      % update iteration number from line search
        update_term = u_IC + ( step_size .* dir_cur );                                                  % retraction operator term
        retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );                   % retraction operator                                   
        u_IC = retraction .* update_term;                                                               % current initial forward state
        manifold_size = sum( u_IC .* conj(u_IC) )*(L1*L2)/N^2;                                          % current manifold (L^2 inner product of initial forward state) 
        IC = 'optimized';                                                                               % change initial condition
        [ v_TC , u_TC ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,u_IC,0);                   % terminal forward state via forward equation
        J_old = J_cur;                                                                                  % old objective functional
        J_cur = sum( u_TC .* conj(u_TC) )*(L1*L2)/N^2;                                                  % current objective functional 
        iter = iter + 1;                                                                                % update iteration number
        projGradJ_old = projGradJ_cur;                                                                  % old projected objective gradient
        dir_old = dir_cur;                                                                              % old direction

        J_history(iter,1) = J_cur;                                                                      % store objective functional history
        J_change(iter,1) = (J_cur - J_old)/(J_old);                                                     % store change in objective functional value
        stepsize_history(iter,1) = step_size;                                                           % store optimization step size
        manifold_history(iter,1) = manifold_size;                                                       % store manifold sizes

    end

    J_change(1,1) = 0;                                                                   % fix initial change in objective functional value 
    %{
    plot_2DKS(save_each, 'kappa', IC, N, dt, T, L_s1, L_s2, kappa,pertIC);                  % save/inspect kappa test figure
    kappalist(:,testcounter) = kappa;                                                       % save kappa test values
    close all                                                                               % close any open figures
    kappalist_file = [pwd '/data/kappa/kappalist_' IC '_p' pertIC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(kappalist, kappalist_file,'Delimiter','tab');
    rieszlist_file = [pwd '/data/kappa/rieszlist_' IC '_p' pertIC '_N_' num2str(N) '' ...
     '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(rieszlist, rieszlist_file,'Delimiter','tab');
    %}
end
function [J_cur , J_history , u_TC , u_IC] = optimize_2DKS(IC,N,L_s1,L_s2,dt,T,u_TC,v_TC,u_IC,Ntime_save_max)

    maxiter = 1000;                                 % set maximum number of iterations
    if exist('J_history','var') == 0
        J_history = NaN(maxiter,1);                        % store updated objective functionals  
        J_change = NaN(maxiter,1);                         % store changes in objective functional value
        stepsize_history = NaN(maxiter,1);                 % store optimization step sizes
        manifold_history = NaN(maxiter,1);                 % store manifold sizes
        time_history = NaN(maxiter,1);                     % store timer count
        gradJsize_history = NaN(maxiter,1);                % store gradient sizes
        momentumsize_history = NaN(maxiter,1);             % store momentum sizes
        diagnostics_history = NaN(maxiter,7);              % store all of the above
        linesearchJ_history = NaN(maxiter,maxiter);        % store line search objective functionals
    end

    originalIC = IC;                                                                        % save name of initial guess IC
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
    gradJsize_history(iter,1) = 0;                                                          % initialize gradient sizes
    time_history(iter,1) = toc;                                                             % initialize timer count
    momentumsize_history(iter,1) = momentum_size;                                           % initialize momentum sizes

    while (abs(J_change(iter,1)) > tol) && (iter <= maxiter)

        disp(['Solving adjoint problem for iteration ' num2str(iter)])
        [~, GradJ] = solve_2DKS(IC,'backward',N,L_s1,L_s2,dt,T,save_each,Ntime_save_max,v_TC,originalIC);        % current objective gradient via adjoint equation
        toc
        GradJ_size = sum( GradJ .* conj(GradJ) )*(L1*L2)/N^2;                                           % current objective gradient size
        angleGradJ = sum( u_IC .* conj(GradJ) )*(L1*L2)/N^2;                                            % angle with current objective gradient
        projGradJ_cur = GradJ - (angleGradJ/manifold_size).*(u_IC);                                     % current projected objective gradient
        if iter > 1
            angleDir_old = sum( u_IC .* conj(dir_old) )*(L1*L2)/N^2;                                    % angle with old direction
            angleprojGradJ_old = sum( u_IC .* conj(projGradJ_old) )*(L1*L2)/N^2;                        % angle with old projected gradient
            vectransport = (dir_old - (angleDir_old/manifold_size).*(u_IC))/manifold_size;              % current vector transport operator
            transportprojGradJ_old = (projGradJ_old - (angleprojGradJ_old/manifold_size).*(u_IC))/manifold_size; % transport old projected gradient 
            diff_projGradJ = projGradJ_cur - transportprojGradJ_old;                                    % momentum parameter term
            diff_projGradJ_size = sum( projGradJ_cur .* conj(diff_projGradJ) )*(L1*L2)/N^2;             % momentum parameter numerator
            projGradJ_old_size = sum( projGradJ_old .* conj(projGradJ_old) )*(L1*L2)/N^2;               % momentum parameter denominator
            momentum_size = diff_projGradJ_size / projGradJ_old_size;                                   % current momentum parameter
        elseif iter == 1
            IC = 'optimized';                                                                           % set IC to optimized
            J_change(1,1) = NaN;                                                                        % fix initial change in objective functional value 
            step_size = 1e-7;%(angleGradJ/GradJ_size);                                                  % initialize current step-size
            stepsize_history(iter,1) = step_size;                                                       % store optimization step size
            diagnostics_history(1,:) = [J_history(1,1), J_change(1,1), stepsize_history(1,1)...
                manifold_history(1,1), time_history(1,1), gradJsize_history(1,1),...
                momentumsize_history(1,1)];
        end
        dir_cur = projGradJ_cur + (momentum_size .* vectransport);                                      % current direction
        [step_size,iter_search,J_search] = optimize_stepsize(dir_cur,u_IC,step_size,IC,N,L_s1,L_s2,dt,T,Ntime_save_max,originalIC); % current step-size via Brent's method
        disp(['Line-search problem number of iterations: ' num2str(iter_search)])
        linesearchJ_history(1:iter_search,iter) = J_search;                                             % store objective functional history from line search
        update_term = u_IC + ( step_size .* dir_cur );                                                  % retraction operator term
        retraction =  sqrt(manifold_history(1,1)) / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );                   % retraction operator                                   
        u_IC = retraction .* update_term;                                                               % current initial forward state
        manifold_size = sum( u_IC .* conj(u_IC) )*(L1*L2)/N^2;                                          % current manifold (L^2 inner product of initial forward state) 
        [ v_TC , u_TC ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,Ntime_save_max,u_IC,originalIC);    % terminal forward state via forward equation
        J_old = J_cur;                                                                                  % old objective functional
        J_cur = sum( u_TC .* conj(u_TC) )*(L1*L2)/N^2;                                                  % current objective functional 
        iter = iter + 1;                                                                                % update iteration number
        projGradJ_old = projGradJ_cur;                                                                  % old projected objective gradient
        dir_old = dir_cur;                                                                              % old direction

        J_history(iter,1) = J_cur;                                                                      % store objective functional history
        J_change(iter,1) = (J_cur - J_old)/(J_old);                                                     % store change in objective functional value
        stepsize_history(iter,1) = step_size;                                                           % store optimization step size
        manifold_history(iter,1) = manifold_size;                                                       % store manifold sizes
        time_history(iter,1) = toc;                                                                     % store timer count 
        gradJsize_history(iter,1) = GradJ_size;                                                         % store gradient sizes
        momentumsize_history(iter,1) = momentum_size;                                                   % store momentum sizes
        diagnostics_history(iter,:) = [J_history(iter,1), J_change(iter,1), stepsize_history(iter,1)...
            manifold_history(iter,1), time_history(iter,1), gradJsize_history(iter,1),...
            momentumsize_history(iter,1)];
    end

    J_history = rmmissing(J_history);
    diagnostics_history = diagnostics_history(1:length(find(~isnan(diagnostics_history(:,1)))),:);      % remove NaN rows
    linesearchJ_history = linesearchJ_history(:,1:(iter-1));                                            % remove NaN columns
    mkdir([pwd  '/data/optimization' ]);
    diagnostics_file = [pwd '/data/optimization/diagnostics_' IC '_' originalIC '_N_' num2str(N) '' ...
    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(diagnostics_history, diagnostics_file,'Delimiter','tab');
    linesearchJ_file = [pwd '/data/optimization/linesearchJ_' IC '_' originalIC '_N_' num2str(N) '' ...
     '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
    writematrix(linesearchJ_history, linesearchJ_file,'Delimiter','tab');
end
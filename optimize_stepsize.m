function [result,iter_search,J_search] = optimize_stepsize(dir_cur,u_IC_cur,step_cur,IC,N,L_s1,L_s2,dt,T)

% find optimal stepsize for objective functional using Brent's method %
    
    iter_search = 0;                                                                           % initialize count of line search iterations
    J_search = NaN(1000,1);                                                                    % initialize storage for objective functional history 
    L1 = 2*pi*L_s1;                                                                            % dimension 1 length
    L2 = 2*pi*L_s2;                                                                            % dimension 2 length
    save_each = T/dt;                                                                          % save only final timestep for forward solver

    TOL = 10^(-5);
    GOLD = 1.618034; GLIMIT = 100;
    ITMAX = 300; CGOLD = .381966; ZEPS = 10^(-10);
    
    % initialize 3 x-points (left, center, right) and keep switching them to converge on minimum value
    
    % This part of the code brackets the location of the minimum of the functional J
    iter = 0;
    tA = 0;                                                                                    % initial left x
    tB = step_cur + 10^(-9);                                                                   % initial center x
    
    update_term = u_IC_cur + ( tA .* dir_cur );                                                % retraction operator term
    retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
    update_term = retraction .* update_term;                                                   % current initial forward state
    [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
    FA = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % initial left f(x) is negative of current objective functional 
    iter_search = iter_search + 1;                                                             % update iteration number
    J_search(iter_search) = -FA;                                                               % store objective functional history
    
    update_term = u_IC_cur + ( tB .* dir_cur );                                                % retraction operator term
    retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
    update_term = retraction .* update_term;                                                   % current initial forward state
    [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
    FB = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % initial center f(x) is negative of current objective functional 
    iter_search = iter_search + 1;                                                             % update iteration number
    J_search(iter_search) = -FB;                                                               % store objective functional history
    
    if FB>FA %% existing code leads to negative tC value
        aux = tA;
        tA = tB; % new center x
        tB = aux; % new right x
        aux = FA;
        FA = FB; % new center f(x)
        FB = aux; % new right f(x)
    end
    
    % FA > FB, tB ? tA
    
    tC = tB + GOLD*(tB - tA);                                                                  % initial right or center x 
    update_term = u_IC_cur + ( tC .* dir_cur );                                                % retraction operator term
    retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
    update_term = retraction .* update_term;                                                   % current initial forward state
    [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
    FC = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % initial right or center f(x) is negative of current objective functional 
    iter_search = iter_search + 1;                                                             % update iteration number
    J_search(iter_search) = -FC;                                                               % store objective functional history
    
    while FB>=FC && iter<ITMAX 
        
        %     SA = (tC-tB)*FA;
        %     SB = (tA-tC)*FB;
        %     SC = (tB-tA)*FC;
        %     tP = 0.5*( (tC+tB)*SA + (tA+tC)*SB + (tB+tA)*SC)/(SA+SB+SC);
        iter = iter+1;
        R = (tB-tA)*(FB-FC); 
        Q = (tB-tC)*(FB-FA); 
        tP = tB - 0.5*((tB-tC)*Q - (tB-tA)*R)/(sign(Q-R)*abs(Q-R));                            % new x correction
        
        Pmax = tB + GLIMIT*(tC-tB); 
    
        if (tB-tP)*(tP-tC)>0                                                                   % if tP is center x
            update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
            retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
            update_term = retraction .* update_term;                                                   % current initial forward state
            [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
            FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
            iter_search = iter_search + 1;                                                             % update iteration number
            J_search(iter_search) = -FP;                                                               % store objective functional history
            
            if FP<FC
                tA = tB;                                                                               % new end x
                FA = FB;                                                                               % new f(x)
                tB = tP;                                                                               % new center x
                FB = FP;                                                                               % new f(x)
                break;
            elseif FP>FB
                tC = tP;
                FC = FP;
                break;
            end
            
            tP = tC + GOLD*(tC-tB);                                                                    % new x correction
            update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
            retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
            update_term = retraction .* update_term;                                                   % current initial forward state
            [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
            FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
            iter_search = iter_search + 1;                                                             % update iteration number
            J_search(iter_search) = -FP;                                                               % store objective functional history
            
        elseif (tC-tP)*(tP-Pmax)>0
            update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
            retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
            update_term = retraction .* update_term;                                                   % current initial forward state
            [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
            FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
            iter_search = iter_search + 1;                                                             % update iteration number
            J_search(iter_search) = -FP;                                                               % store objective functional history
            
            if FP<FC
                tB = tC;
                tC = tP;
                FB = FC;
                FC = FP;
                tP = tC+GOLD*(tC-tB);                                                                      % new x correction
                update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
                retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
                update_term = retraction .* update_term;                                                   % current initial forward state
                [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
                FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
                iter_search = iter_search + 1;                                                             % update iteration number
                J_search(iter_search) = -FP;                                                               % store objective functional history
            end
            
        elseif (tP-Pmax)*(Pmax-tC)>=0
            tP = Pmax;                                                                                 % new x correction
            update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
            retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
            update_term = retraction .* update_term;                                                   % current initial forward state
            [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
            FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
            iter_search = iter_search + 1;                                                             % update iteration number
            J_search(iter_search) = -FP;                                                               % store objective functional history        
        else
            tP = tC + GOLD*(tC-tB);                                                                    % new x correction
            update_term = u_IC_cur + ( tP .* dir_cur );                                                % retraction operator term
            retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
            update_term = retraction .* update_term;                                                   % current initial forward state
            [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
            FP = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
            iter_search = iter_search + 1;                                                             % update iteration number
            J_search(iter_search) = -FP;                                                               % store objective functional history
        end
        
        tA = tB;
        tB = tC;
        tC = tP;
        FA = FB;
        FB = FC;
        FC = FP;
    
    end
    
    if iter==ITMAX
        result = step_cur;
        return
    end
    
    % This part of the code finds the minimum of the functional and gives the value of the optimal step size

    D = 0;
    A = min(tA,tC);
    B = max(tA,tC);
    V = tB; W = V; X = V; E = 0;
    FX = FB; FV = FX; FW = FX;
  
    for j=1:ITMAX
              
        XM = 0.5*(A+B);                                                                                % halfway between [A,B]
        TOL1 = (TOL*abs(X)+ZEPS);
        TOL2 = 2*TOL1;
        x_step = abs(X-XM);
        x_tol = (TOL2-0.5*(B-A));
    
        if ( x_step <= x_tol )
            break;
        end

        FLAG = 1;
        if ( abs(E) > TOL1 )
            R = (X-W)*(FX-FV);
            Q = (X-V)*(FX-FW);
            P = (X-V)*Q - (X-W)*R;
            Q = 2*(Q-R);
            if Q>0
                P = -P;
            end
            Q = abs(Q);
            ETEMP = E;
            E = D;
            
            if (abs(P) >= abs(0.5*Q*ETEMP)) || (P <= Q*(A-X)) || ...
                    (P >= Q*(B-X))
                FLAG = 1;
            else
                FLAG = 2;
            end
            
        end
        
        switch FLAG
            case 1
                if X >= XM
                    E = A-X;
                else
                    E=B-X;
                end
                D = CGOLD*E;
            case 2
                D = P/Q;
                U = X+D;
                if (U-A < TOL2) || (B-U < TOL2)
                    D = sign(XM-X)*TOL1;
                end
        end
        
        if abs(D) >= TOL1
            U = X+D;
        else
            U = X + sign(D)*TOL1;
        end
        
        update_term = u_IC_cur + ( U .* dir_cur );                                                 % retraction operator term
        retraction =  1 / sqrt(sum( update_term .* conj(update_term) )*(L1*L2)/N^2 );              % retraction operator                                   
        update_term = retraction .* update_term;                                                   % current initial forward state
        [ ~ , u_TC_cur ] = solve_2DKS(IC,'forward',N,L_s1,L_s2,dt,T,save_each,update_term,0);      % terminal forward state via forward equation
        FU = - ( sum( u_TC_cur .* conj(u_TC_cur) )*(L1*L2)/N^2 );                                  % new f(x) is negative of current objective functional 
        iter_search = iter_search + 1;                                                             % update iteration number
        J_search(iter_search) = -FU;                                                               % store objective functional history
        
        if FU <= FX
            if U >= X
                A = X;
            else
                B = X;
            end
            V = W;
            FV = FW;
            W = X;
            FW = FX;
            X = U;
            FX = FU;
        else
            if U < X
                A = U;
            else
                B = U;
            end
            
            if (FU <= FW) || (W == X)
                V = W;
                FV = FW;
                W = U;
                FW = FU;
            elseif (FU <= FV) || (V==X) || (V==W)
                V = U;
                FV = FU;
            end
        end
    
    end

    result = X;
    J_search = rmmissing(J_search);

return
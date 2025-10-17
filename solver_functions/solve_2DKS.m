function [ v_out , u_out , utilityout ] = solve_2DKS(IC,solver,N,K,L_s1,L_s2,dt,T,save_each,Ntime_save_max,utility1,utility2)

%{ 
Output:
v_n: solution vector for each time step in Fourier space
u_n: solution vector for each time step in Physical space
%}

%%% (1) initialize problem %%%

    N_x2 = N;                                                % discretized equally in each dimension

    % length-scale parameters
    L_x1 = (1/L_s1);
    L_x2 = (1/L_s2);
    L1 = 2*pi*L_s1;
    L2 = 2*pi*L_s2;

    % unit physical space domain
    x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
    x2_pts = L2*linspace( 0 , 1 - 1/N_x2 , N_x2 ); 
    [ x1 , x2 ] = meshgrid(x1_pts,x2_pts);                      % 2-dimensional grid
    
    % fourier space domain for nonlinear term
    k1_n_pts = 2*pi/L1*[ 0 : N/2-1 , 0 , -N/2+1 : -1]; 
    k2_n_pts = 2*pi/L2*[ 0 : N_x2/2-1 , 0 , -N_x2/2+1 : -1]; 
    [ k1_n , k2_n ] = meshgrid(k1_n_pts,k2_n_pts);              % 2-dimensional grid

    % fourier space domain for linear term
    k1_l_pts = 2*pi/L1*[0:N/2 -N/2+1:-1];
    k2_l_pts = 2*pi/L2*[0:N_x2/2 -N_x2/2+1:-1];
    [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts);              % 2-dimensional grid

    % reshape
    k1vecN = k1_n(:);
    k2vecN = k2_n(:);
    Klap = k1_l.^2 + k2_l.^2;                                   % for 2nd order linear term (Laplacian) 
    kvecl2 = Klap(:);
    KK = (k1_l.^2 + k2_l.^2).^2;                                % for 4th order linear term (bi-Laplacian)
    kvecl4 = KK(:);                           
    
    % Fourier space operators
    Lin = 1i^2 * (kvecl2) + 1i^4 * (kvecl4);                    % Linear operator
    Lap = 1i^2 * (kvecl2);                                      % Laplace operator
    D1vec = 1i*k1vecN;                                          % Differential operator
    D2vec = 1i*k2vecN;                                          % Differential operator
    
    % coefficient terms (Alimo et al. 2021, Table 2)
    alpha_I = [ 343038331393/1130875731271 , 288176579239/1140253497719 , 253330171251/677500478386 , 189462239225/1091147436423 ];
    beta_I = [ 35965327958/140127563663 , 19632212512/2700543775099 , -173747147147/351772688865 , 91958533623/727726057489 ];
    alpha_E = [ 14/25 , 777974228744/1346157007247 , 251277807242/1103637129625 , 113091689455/220187950967 ];
    beta_E = [ 0 , -251352885992/790610919619 , -383714262797/1103637129625 , -403360439203/1888264787188 ];
    
    % impose initial and boundary conditions in physical and fourier space
    switch IC 
        case 'optimized'
            u_0 = reshape( utility1, [ N , N_x2 ] );
        case 's'
            u_0 = sin( (x1 + x2) ) + sin( x1 ) + sin( x2 );
        case 's1'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) + sin( L_x2*x2 );
        case 'sc1'
            u_0 = cos( L_x1*x1 ) .* sin( L_x2*x2 );
        case 's30'
            u_0 = sin( 1*(L_x1*x1 + L_x2*x2) ) + sin( 30*L_x1*x1 ) + sin( 30*L_x2*x2 );
        case 'tg1'
            u_0 = sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'tg30'
            u_0 = sin( L_x1*x1 ) .* sin( 30*L_x2*x2 );
        case 'stg1'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'stg30'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) +  sin( L_x1*x1 ) .* sin( 30*L_x2*x2 );
        case 'gauss'
            u_0 = exp(-0.1*( (x1 - 0.5*L1).^2 + (x2 - 0.5*L2).^2));
        case 'noise'
            u_0 = 1e-14*(2*rand(N)-ones(N));
        case 'noise4'
            u_0 = 1e-4*(2*rand(N)-ones(N));
        case 'mn5'
            u_0 = 1e-5*( sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                                 sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                                 sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ) ;
    end

    % set initial L^2 energy magnitude
    u_0_mag = sqrt(sum( u_0(:) .* conj(u_0(:)) )*(L1*L2)/N^2);                  % compute norm of IC
    if K ~= 0 
        u_0 = sqrt(K)*(u_0/u_0_mag);                                                    % set norm of IC equal to K
    end
    
    % perturbation function for adjoint calculus
    switch utility2 
        case 's'
            u_pert = sin( (x1 + x2) ) + sin( x1 ) + sin( x2 );
        case 's1'
            u_pert = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) + sin( L_x2*x2 );
        case 's30'
            u_pert = sin( 1*(L_x1*x1 + L_x2*x2) ) + sin( 30*L_x1*x1 ) + sin( 30*L_x2*x2 );
        case 'tg1'
            u_pert = sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'tg30'
            u_pert = sin( L_x1*x1 ) .* sin( 30*L_x2*x2 );
        case 'stg1'
            u_pert = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'stg30'
            u_pert = sin( (L_x1*x1 + L_x2*x2) ) +  sin( L_x1*x1 ) .* sin( 30*L_x2*x2 );
    end                   

    switch solver
        case 'forward'
            v_0 = fft2(u_0);            % FFT of physical initial condition
            v_out = 0;
            u_out = 0;
            utilityout = u_0(:);
        case 'kappa'
            eps = utility1;            % perturbation magnitude for kappa test
            u_0 = u_0 + eps*u_pert;     % perturbed initial condition
            v_0 = fft2(u_0);            % FFT of perturbed physical initial condition
            v_out = 0;
            u_out = 0;
            utilityout = 0;
        case 'backward'
            try
                utilityout = u_pert(:);     % perturbed physical initial condition
            catch
                utilityout = 0;
            end
            v_TC = utility1;           % forward solution in Fourier space
    end

%%% (2) solve time-dependent problem %%%

    % number of timesteps
    time1 = ceil(T/dt); 
    time2 = ceil(time1/save_each);
    if T >= 0
        Ntime = max(time1,time2);
        Ntime_save = min(time1,time2);
        save_each = ceil(Ntime/Ntime_save);
        Ntime_save = ceil(Ntime/save_each);
        Ntime = save_each*Ntime_save;
    else
        Ntime_save = 1;
        Ntime = time1;
        save_each = time1;
    end

    % solution vectors in fourier and physical spectrum (max 10k columns)
    switch solver
        case {'forward','kappa'}  
            if Ntime_save > Ntime_save_max
                u_n = zeros( N * N_x2 , Ntime_save_max );
                v_n = zeros( N * N_x2 , Ntime_save_max );
            else
                u_n = zeros( N * N_x2 , Ntime_save );
                v_n = zeros( N * N_x2 , Ntime_save );
            end
        case {'backward'}
            u_n = zeros( N * N_x2 , 2 );
            v_n = zeros( N * N_x2 , 2 );
    end

    switch solver 
        case {'forward'}                                                                       % Solve vectorized equation by IMEXRK4 method
            u_n(:,1) = u_0(:);                                                                  % physical IC
            v_n(:,1) = v_0(:);                                                                  % fourier IC
            v_step = v_n(:,1);                                                                  % initialize stepping with fourier IC
            Nonlin_v0 = 0;    
            for i = 2:Ntime

                % nonlinear terms and solution substeps
                for k = 1:4
                    v_step2x = reshape( D1vec .* v_step, [ N , N_x2 ] );                     % f_x
                    v_step2y = reshape( D2vec .* v_step, [ N , N_x2 ] );                     % f_y
                    w1_r = multiply2D( v_step2x , v_step2x , 'fourier2real' );                  % f_x * f_x in physical space (pseudospectral)
                    w1s_r = multiply2D( v_step2y , v_step2y , 'fourier2real' );                 % f_y * f_y in physical space (pseudospectral)
                    Nonlin_v1_r = (1/2) * ( w1_r + w1s_r );                                     % (1/2)*(f_x * f_x + f_y * f_y) in physical space 
                    Nonlin_v1 = multiply2D( fft2(Nonlin_v1_r) , fft2(Nonlin_v1_r) ,'dealias');  % dealias
                    Nonlin_v1 = Nonlin_v1(:);
    
                    v_1 = ( 1 + (dt * alpha_I(k) * Lin) ).^(-1) .* ...
                        ( ( 1 - (dt * beta_I(k) * Lin) ) .* v_step - ...
                        (dt * alpha_E(k) * Nonlin_v1) - (dt * beta_E(k) * Nonlin_v0) );

                    v_step = v_1;
                    Nonlin_v0 = Nonlin_v1;
                end

                % correction of non-zero mean solution
                if abs(mean(v_1)) > 1e-5 
                    v_1(1) = 0;
                end

                % save solution step to workspace
                v_12 = reshape( v_1, [ N , N_x2 ] );
                u_1 = real(ifft2(v_12));  
                v_step = v_1;
                
                if (save_each == 1) && (mod(i,Ntime_save_max) ~= 0)
                    v_n(:,mod(i,Ntime_save_max)) = v_12(:);
                    u_n(:,mod(i,Ntime_save_max)) = u_1(:);
                elseif (i == Ntime) || ( (save_each == 1) && (mod(i,Ntime_save_max) == 0) )
                    v_n(:,end) = v_12(:);
                    u_n(:,end) = u_1(:);
                elseif mod(i,save_each) == 0
                    v_n(:,mod(i/save_each + 1,Ntime_save_max)) = v_12(:);
                    u_n(:,mod(i/save_each + 1,Ntime_save_max)) = u_1(:);
                end

                if (Ntime_save > Ntime_save_max) && (mod(i,Ntime_save_max) == 0) && (i ~= Ntime)
                    currentT = i/Ntime*T;
                    time = toc;
                    save_2DKSsolution('forward', u_n, v_n, time, IC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, utility2); % save solution to machine
                    if (Ntime - i) < Ntime_save_max
                        fullsave = 0;
                        Ntime_save_remaining = Ntime - i;
                        u_n = zeros( N * N_x2 , Ntime_save_remaining );
                        v_n = zeros( N * N_x2 , Ntime_save_remaining );
                    else
                        fullsave = 1;
                        u_n = zeros( N * N_x2 , Ntime_save_max );
                        v_n = zeros( N * N_x2 , Ntime_save_max );
                    end
                    v_out = v_n(:,end);
                    u_out = u_n(:,end);
                elseif (Ntime_save > Ntime_save_max) && (i == Ntime)
                    time = toc;
                    if fullsave == 1
                        save_2DKSsolution('forward', u_n, v_n, time, IC, dt, T, N, K, L_s1, L_s2, Ntime_save_max, utility2); % save solution to machine
                    else 
                        save_2DKSsolution('forward', u_n, v_n, time, IC, dt, T, N, K, L_s1, L_s2, Ntime_save_remaining, utility2); % save solution to machine
                    end
                    v_out = v_n(:,end);
                    u_out = u_n(:,end);
                elseif i == Ntime % only reached for single-file evolution data
                    time = toc;
                    save_2DKSsolution('forward', u_n, v_n, time, IC, dt, T, N, K, L_s1, L_s2, Ntime_save, utility2); % save solution to machine
                    v_out = v_n(:,end);
                    u_out = u_n(:,end);
                end
            end
        case {'kappa'}                                                                % Solve vectorized equation by IMEXRK4 method
            u_n(:,1) = u_0(:);                                                                  % physical IC
            v_n(:,1) = v_0(:);                                                                  % fourier IC
            v_step = v_n(:,1);                                                                  % initialize stepping with fourier IC
            Nonlin_v0 = 0;    
            for i = 2:Ntime

                % nonlinear terms and solution substeps
                for k = 1:4
                    v_step2x = reshape( D1vec .* v_step, [ N , N_x2 ] );                     % f_x
                    v_step2y = reshape( D2vec .* v_step, [ N , N_x2 ] );                     % f_y
                    w1_r = multiply2D( v_step2x , v_step2x , 'fourier2real' );                  % f_x * f_x in physical space (pseudospectral)
                    w1s_r = multiply2D( v_step2y , v_step2y , 'fourier2real' );                 % f_y * f_y in physical space (pseudospectral)
                    Nonlin_v1_r = (1/2) * ( w1_r + w1s_r );                                     % (1/2)*(f_x * f_x + f_y * f_y) in physical space 
                    Nonlin_v1 = multiply2D( fft2(Nonlin_v1_r) , fft2(Nonlin_v1_r) ,'dealias');  % dealias
                    Nonlin_v1 = Nonlin_v1(:);
    
                    v_1 = ( 1 + (dt * alpha_I(k) * Lin) ).^(-1) .* ...
                        ( ( 1 - (dt * beta_I(k) * Lin) ) .* v_step - ...
                        (dt * alpha_E(k) * Nonlin_v1) - (dt * beta_E(k) * Nonlin_v0) );

                    v_step = v_1;
                    Nonlin_v0 = Nonlin_v1;
                end

                % correction of non-zero mean solution
                if abs(mean(v_1)) > 1e-5 
                    v_1(1) = 0;
                end

                % save solution step to workspace
                v_12 = reshape( v_1, [ N , N_x2 ] );
                u_1 = real(ifft2(v_12));
                v_step = v_1;

                if i == Ntime
                    v_n(:,end) = v_12(:);
                    u_n(:,end) = u_1(:);
                elseif mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_12(:);
                    u_n(:,i/save_each + 1) = u_1(:);
                end

                if i == Ntime
                    v_out = v_n(:,end);
                    u_out = u_n(:,end);
                end
            end
        case {'backward'} % Solve adjoint equation
            v_n(:,2) = 2 .* v_TC;                                                           % fourier TC
            u_n(:,2) = real(ifft2(v_n(:,2)));          
            v_step = v_n(:,end);                                                             % initialize stepping with fourier TC
            Nonlin_v0 = 0;
            Ntime_residual = mod(Ntime,Ntime_save_max);
            if Ntime_residual == 0
                Ntime_residual = Ntime_save_max;
            end
            [~, v_fwd] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_residual, utility2);
            count = 0;
            for i = Ntime-1:-1:1
                v_fwdstep = v_fwd(:,Ntime_residual-count);
                v_fwd2x = reshape( D1vec .* v_fwdstep, [ N , N_x2 ] );                       % f_x
                v_fwd2y = reshape( D2vec .* v_fwdstep, [ N , N_x2 ] );                       % f_y
                v_fwd2lap = reshape( Lap .* v_fwdstep, [ N , N_x2 ] );                       % lap(f)

                % nonlinear terms and solution substeps
                for k = 1:4
                    v_step2 = reshape( v_step, [ N , N_x2 ] );                               % z (adjoint variable)
                    v_step2x = reshape( D1vec .* v_step, [ N , N_x2 ] );                     % z_x
                    v_step2y = reshape( D2vec .* v_step, [ N , N_x2 ] );                     % z_y
                    w1_r = multiply2D( v_step2x , v_fwd2x , 'fourier2real' );                   % z_x * f_x in physical space (pseudospectral)
                    w1s_r = multiply2D( v_step2y , v_fwd2y , 'fourier2real' );                  % z_y * f_y in physical space (pseudospectral)
                    w2_r = multiply2D( v_fwd2lap , v_step2 , 'fourier2real' );                  % lap(f) * z in physical space (pseudospectral)
                    Nonlin_v1_r = (-1)*( w1_r + w1s_r ) - ( w2_r ) ;                            % - (z_x * f_x + z_y * f_y) - (lap(f) * z) in physical space
                    Nonlin_v1 = multiply2D( fft2(Nonlin_v1_r) , fft2(Nonlin_v1_r) , 'dealias'); % dealias
                    Nonlin_v1 = Nonlin_v1(:);

                    v_1 = ( 1 + (dt * alpha_I(k) * Lin) ).^(-1) .* ...
                        ( ( 1 - (dt * beta_I(k) * Lin) ) .* v_step - ...
                        (dt * alpha_E(k) * Nonlin_v1) - (dt * beta_E(k) * Nonlin_v0) );

                    v_step = v_1;
                    Nonlin_v0 = Nonlin_v1;
                end

                % correction of non-zero mean solution
                if abs(mean(v_1)) > 1e-5 
                    v_1(1) = 0;
                end

                % save solution step to workspace
                v_12 = reshape( v_1, [ N , N_x2 ] );
                u_1 = real(ifft2(v_12));
                v_step = v_1;

                if i == 1
                    time = toc;
                    v_n(:,1) = v_step;
                    u_n(:,1) = real(ifft2(v_n(:,1)));
                    save_2DKSsolution('backward', u_n, v_n, time, IC, dt, T, N, K, L_s1, L_s2, 2, utility2); % save solution to machine
                    v_out = v_12(:);
                    u_out = u_1(:);
                elseif mod(count,Ntime_residual) == (Ntime_residual - 1)
                    Ntime_residual = Ntime_save_max;
                    count = 0;
                    currentT = i/Ntime*T;
                    [~, v_fwd] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, Ntime_residual, utility2);
                else
                    count = count + 1;                    
                end
            end

    end
end
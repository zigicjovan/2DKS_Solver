function [ v_n , u_n , utility ] = solve_2DKS(IC,solver,N_x1,L_s1,L_s2,dt,T,save_each,utility)

%{ 
Output:
v_n: solution vector for each time step in Fourier space
u_n: solution vector for each time step in Physical space
%}

%%% (1) initialize problem %%%

    N_x2 = N_x1;                                                % discretized equally in each dimension

    % length-scale parameters
    L_x1 = (1/L_s1);
    L_x2 = (1/L_s2);
    L1 = 2*pi*L_s1;
    L2 = 2*pi*L_s2;

    % unit physical space domain
    x1_pts = L1*linspace( 0 , 1 - 1/N_x1 , N_x1 ); 
    x2_pts = L2*linspace( 0 , 1 - 1/N_x2 , N_x2 ); 
    [ x1 , x2 ] = meshgrid(x1_pts,x2_pts);                      % 2-dimensional grid
    
    % fourier space domain for nonlinear term
    k1_n_pts = 2*pi/L1*[ 0 : N_x1/2-1 , 0 , -N_x1/2+1 : -1]; 
    k2_n_pts = 2*pi/L2*[ 0 : N_x2/2-1 , 0 , -N_x2/2+1 : -1]; 
    [ k1_n , k2_n ] = meshgrid(k1_n_pts,k2_n_pts);              % 2-dimensional grid

    % fourier space domain for linear term
    k1_l_pts = 2*pi/L1*[0:N_x1/2 -N_x1/2+1:-1];
    k2_l_pts = 2*pi/L2*[0:N_x2/2 -N_x2/2+1:-1];
    [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts);              % 2-dimensional grid

    % reshape
    k1vecN = k1_n(:);
    k2vecN = k2_n(:);
    K = k1_l.^2 + k2_l.^2;                                      % for 2nd order linear term (Laplacian) 
    kvecl2 = K(:);
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
        case 'sin'
            u_0 = sin( (x1 + x2) ) + sin( x1 ) + sin( x2 );
        case 'sinL'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) + sin( L_x2*x2 );
        case 'gauss'
            u_0 = exp(-0.1*( (x1 - 0.5*L1).^2 + (x2 - 0.5*L2).^2));
        case 'noise'
            u_0 = 1e-14*(2*rand(N_x1)-ones(N_x1));
        case 'mn1'
            u_0 = ((10.0)^(-1))*( sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                                 sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                                 sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ) ;
        case 'mn5'
            u_0 = ((10.0)^(-5))*( sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                                 sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                                 sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ) ;
        case 'mn10'
            u_0 = ((10.0)^(-10))*( sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                                  sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                                  sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ) ;
        case 'mn14'
            u_0 = ((10.0)^(-14))*( sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                                  sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                                  sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ) ;
    end
    switch solver
        case 'forward'
            v_0 = fft2(u_0); 
        case 'kappa'
            eps = utility;
            u_pert = u_0;
            u_0 = u_0 + eps*u_pert;
            v_0 = fft2(u_0); 
        case 'backward'
            v_fwd = utility;
    end

%%% (2) solve time-dependent problem %%%

    % number of timesteps
    time1 = ceil(T/dt); 
    time2 = ceil(time1/save_each);
    if T >= 1
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

    % solution vectors in fourier and physical spectrum
    u_n = zeros( N_x1 * N_x2 , Ntime_save );
    v_n = zeros( N_x1 * N_x2 , Ntime_save );

    switch solver 
        case {'forward','kappa'}                                                                % Solve vectorized equation by IMEXRK4 method
            u_n(:,1) = u_0(:);                                                                  % physical IC
            v_n(:,1) = v_0(:);                                                                  % fourier IC
            v_step = v_n(:,1);                                                                  % initialize stepping with fourier IC
            Nonlin_v0 = 0;    
            for i = 2:Ntime

                % nonlinear terms and solution substeps
                for k = 1:4
                    v_step2x = reshape( D1vec .* v_step, [ N_x1 , N_x2 ] );                     % f_x
                    v_step2y = reshape( D2vec .* v_step, [ N_x1 , N_x2 ] );                     % f_y
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
                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                u_1 = real(ifft2(v_12));  
                v_step = v_1;
                
                if save_each == 1
                    v_n(:,i) = v_12(:);
                    u_n(:,i) = u_1(:);
                elseif mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_12(:);
                    u_n(:,i/save_each + 1) = u_1(:);
                elseif i == Ntime
                    v_n(:,Ntime) = v_12(:);
                    u_n(:,Ntime) = u_1(:);
                end

            end
        case {'backward'} % Solve adjoint equation
            v_n(:,end) = 2 .* v_fwd(:,end);                                                     % fourier TC
            u_n(:,end) = real(ifft2(v_n(:,end)));          
            v_step = v_n(:,end);                                                                % initialize stepping with fourier TC
            Nonlin_v0 = 0;
            for i = Ntime-1:-1:1
                v_fwdstep = v_fwd(:,i);

                % nonlinear terms and solution substeps
                for k = 1:4
                    v_step2 = reshape( v_step, [ N_x1 , N_x2 ] );                               % z (adjoint variable)
                    v_step2x = reshape( D1vec .* v_step, [ N_x1 , N_x2 ] );                     % z_x
                    v_step2y = reshape( D2vec .* v_step, [ N_x1 , N_x2 ] );                     % z_y
                    v_fwd2x = reshape( D1vec .* v_fwdstep, [ N_x1 , N_x2 ] );                   % f_x
                    v_fwd2y = reshape( D2vec .* v_fwdstep, [ N_x1 , N_x2 ] );                   % f_y
                    v_fwd2lap = reshape( Lap .* v_fwdstep, [ N_x1 , N_x2 ] );                   % lap(f)
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
                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                u_1 = real(ifft2(v_12));  
                v_step = v_1;
                
                if save_each == 1
                    v_n(:,i) = v_12(:);
                    u_n(:,i) = u_1(:);
                elseif mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_12(:);
                    u_n(:,i/save_each + 1) = u_1(:);
                elseif i == Ntime
                    v_n(:,Ntime) = v_12(:);
                    u_n(:,Ntime) = u_1(:);
                end

            end

    end
end
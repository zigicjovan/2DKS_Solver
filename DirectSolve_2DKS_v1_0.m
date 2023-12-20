function [ v_n , u_n ] = DirectSolve_2DKS_v1_0(IC,method,N_x1,L_s1,L_s2,dt,T,save_each)

%{ 
Output:
v_n: solution vector for each time step in Fourier space
u_n: solution vector for each time step in Physical space
%}

%%% (1) initialize problem %%%

    % Discretized equally in each dimension
    N_x2 = N_x1;

    % length-scale parameters
    % L_x1 = (pi/L_s1)^2;
    % L_x2 = (pi/L_s2)^2;
    L_x1 = (1/L_s1);
    L_x2 = (1/L_s2);
    L1 = 2*pi/L_x1;
    L2 = 2*pi/L_x2;

    % unit physical space domain
    x1_pts = L1*linspace( 0 , 1 - 1/N_x1 , N_x1 ); 
    x2_pts = L2*linspace( 0 , 1 - 1/N_x2 , N_x2 ); 
    [ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
    % fourier space domain for nonlinear term
    k1_n_pts = [ 0 : N_x1/2-1 , 0 , -N_x1/2+1 : -1]; 
    k2_n_pts = [ 0 : N_x2/2-1 , 0 , -N_x2/2+1 : -1]; 
    [ k1_n , k2_n ] = meshgrid(k1_n_pts,k2_n_pts); % 2-dimensional grid

    % fourier space domain for linear term
    k1_l_pts = [0:N_x1/2 -N_x1/2+1:-1];
    k2_l_pts = [0:N_x2/2 -N_x2/2+1:-1];
    [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts); % 2-dimensional grid

    % reshape
    k1vecN = k1_n(:);
    k2vecN = k2_n(:);
    K = k1_l.^2 + k2_l.^2;
    kvecl2 = K(:);
    % KK = k1_1.^4 + 2*(k1_1.^2 .* k2_1.^2) + k2_1.^4;
    KK = k1_l.^4 + k2_l.^4;
    kvecl4 = KK(:);
    
    switch method
        case 'imexrk4vec'
            % Linear Operator in Fourier Space
            Lin = (-1) * ( 1i^2 * (kvecl2) + 1i^4 * (kvecl4) ); 
        
            % Differential Operators in Fourier Space
            D1vec = 1i*2*pi*k1vecN; 
            D2vec = 1i*2*pi*k2vecN; 
            
            % coefficient terms (Alimo et al. 2021, Table 2)
            alpha_I = [ 343038331393/1130875731271 , 288176579239/1140253497719 ,...
                253330171251/677500478386 , 189462239225/1091147436423 ];
            beta_I = [ 35965327958/140127563663 , 19632212512/2700543775099 ,...
                -173747147147/351772688865 , 91958533623/727726057489 ];
            alpha_E = [ 14/25 , 777974228744/1346157007247 ,...
                251277807242/1103637129625 , 113091689455/220187950967 ];
            beta_E = [ 0 , -251352885992/790610919619 ,...
                -383714262797/1103637129625 , -403360439203/1888264787188 ];
    end
    
    % impose initial and boundary conditions in physical and fourier space
    switch IC 
        case 'sin'
            u_0 = sin( (x1 + x2) ) + sin( x1 ) + sin( x2 );
        case 'gauss'
            u_0 = exp(-0.1*( (x1 - 0.5*L1).^2 + (x2 - 0.5*L2).^2));
    end
    u_0(end,:) = u_0(1,:); % periodic boundary conditions
    u_0(:,end) = u_0(:,1); % periodic boundary conditions
    v_0 = fft2(u_0); 

%%% (2) solve time-dependent problem %%%

    % number of timesteps (clean this up eventually...)
    time1 = ceil(T/dt); 
    time2 = ceil(time1/save_each);
    if T >= 1
        Ntime = max(time1,time2);
        Ntime_save = min(time1,time2);
        save_each = ceil(Ntime/Ntime_save);
        Ntime_save = ceil(Ntime/save_each);
        Ntime = save_each*Ntime_save;
    else
        Ntime_save = time1;
        Ntime = time1;
    end

    % solution vectors in fourier and physical spectrum
    u_n = zeros( N_x1 * N_x2 , Ntime_save );
    u_n(:,1) = u_0(:);
    v_n = zeros( N_x1 * N_x2 , Ntime_save );
    v_n(:,1) = v_0(:);

    %{
    % inspect physical and fourier spectrum
    v_i = reshape( v_n(:,1) , [ N_x1 , N_x2 ] ); 
    u_i = reshape( u_n(:,1) , [ N_x1 , N_x2 ] );
    figure(1); 
    surf(k1_0,k2_0,abs(v_i)); set(gca,'xscale','log','yscale','log','zscale','log');
    % view(90,0);
    figure(2);
    surf(x1,x2,u_i); 
    %}

    switch method % 'imexrk4vec' is the only one guaranteed to work - all others are under construction

        case 'imexrk4vec' % Solve vectorized equation by IMEXRK4 method
            v_0step = v_n(:,1);
            for i = 2:Ntime

                %
                % nonlinear terms and solution substeps
                v_0step2 = reshape( v_0step, [ N_x1 , N_x2 ] );
                w1 = multiply2D( v_0step2 , v_0step2 , 'fourier' );
                Nonlin_v0 = (-1/2) * ( D1vec .* w1(:) + D2vec .* w1(:) );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );

                v_s12 = reshape( v_s1, [ N_x1 , N_x2 ] );
                w2 = multiply2D( v_s12 , v_s12 , 'fourier' );
                Nonlin_v1 = (-1/2) * ( D1vec .* w2(:) + D2vec .* w2(:) );
                v_s2 = ( 1 - (dt * alpha_I(2) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(2) * Lin) ) .* v_s1 + ...
                    (dt * alpha_E(2) * Nonlin_v1) + ...
                    (dt * beta_E(2)) * Nonlin_v0);

                v_s22 = reshape( v_s2, [ N_x1 , N_x2 ] );
                w3 = multiply2D( v_s22 , v_s22 , 'fourier' );
                Nonlin_v2 = (-1/2) * ( D1vec .* w3(:) + D2vec .* w3(:) );
                v_s3 = ( 1 - (dt * alpha_I(3) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(3) * Lin) ) .* v_s2 + ...
                    (dt * alpha_E(3) * Nonlin_v2) + ...
                    (dt * beta_E(3) * Nonlin_v1) ); 

                v_s32 = reshape( v_s3, [ N_x1 , N_x2 ] );
                w4 = multiply2D( v_s32 , v_s32 , 'fourier' );
                Nonlin_v3 = (-1/2) * ( D1vec .* w4(:) + D2vec .* w4(:) );
                v_1 = ( 1 - (dt * alpha_I(4) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(4) * Lin) ) .* v_s3 + ...
                    (dt * alpha_E(4) * Nonlin_v3) + ...
                    (dt * beta_E(4) * Nonlin_v2) ); 
                %}

                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                %{
                figure(1);
                surfc(x1,x2,real(ifft2(v_0step2)))
                figure(2);
                surfc(x1,x2,real(ifft2(v_s12)))
                figure(3);
                surfc(x1,x2,real(ifft2(v_s22)))
                figure(4);
                surfc(x1,x2,real(ifft2(v_s32)))
                figure(5);
                surfc(x1,x2,real(ifft2(v_12)))
                %}
                u_1 = real(ifft2(v_12));  
                u_1(end,:) = u_1(1,:); % boundary conditions
                u_1(:,end) = u_1(:,1); % boundary conditions
                v_0step2 = fft2(u_1);
                v_0step = v_0step2(:);

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
        case 'cnrkw3vec' % Solve vectorized equation by CNRKW3 method
            v_0step = v_n(:,1);
            for i = 2:Ntime
                
                % coefficient terms (Le and Moin 1991)
                alpha_I = [ 4/15 , 1/15 , 1/6 ];
                beta_I = [ 4/15 , 1/15 , 1/6 ];
                alpha_E = [ 8/15 , 5/12 , 3/4 ];
                beta_E = [ 0 , -17/60 , -5/12 ];

                % Differential Operators in Fourier Space
                D1vec = 1i*2*pi*k1vecN; 
                D2vec = 1i*2*pi*k2vecN; 

                % linear terms
                Lin = (-1) * ( 1i^2 * (kvecl2) + 1i^4 * (kvecl4) );  

                %
                % nonlinear terms and solution substeps
                v_0step2 = reshape( v_0step, [ N_x1 , N_x2 ] );
                w = multiply2D( v_0step2 , v_0step2 , 'fourier' );
                Nonlin_v0 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );

                v_s12 = reshape( v_s1, [ N_x1 , N_x2 ] );
                w = multiply2D( v_s12 , v_s12 , 'fourier' );
                Nonlin_v1 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_s2 = ( 1 - (dt * alpha_I(2) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(2) * Lin) ) .* v_s1 + ... 
                    (dt * alpha_E(2) * Nonlin_v1) + ...
                    (dt * beta_E(2) * Nonlin_v0) );

                v_s22 = reshape( v_s2, [ N_x1 , N_x2 ] );
                w = multiply2D( v_s22 , v_s22 , 'fourier' );
                Nonlin_v2 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_1 = ( 1 - (dt * alpha_I(3) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(3) * Lin) ) .* v_s2 + ...
                    (dt * alpha_E(3) * Nonlin_v2) + ...
                    (dt * beta_E(3) * Nonlin_v1) ); 
                %}

                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                %{
                figure(1);
                surfc(x1,x2,real(ifft2(v_0step2)))
                figure(2);
                surfc(x1,x2,real(ifft2(v_s12)))
                figure(3);
                surfc(x1,x2,real(ifft2(v_s22)))
                figure(4);
                surfc(x1,x2,real(ifft2(v_12)))
                %}
                u_1 = real(ifft2(v_12));  
                u_1(end,:) = u_1(1,:); % boundary conditions
                u_1(:,end) = u_1(:,1); % boundary conditions
                v_0step2 = fft2(u_1);
                v_0step = v_0step2(:);

                if save_each == 1
                    v_n(:,i) = v_12(:);
                    u_n(:,i) = u_1(:);
                elseif mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_12(:);
                    u_n(:,i/save_each + 1) = u_1(:);
                end

            end
        case 'imexbdvec' % Solve vectorized equation by IMEXBDF method
            v_0 = v_n(:,1);
            v_1 = v_n(:,1);
            % v12 = reshape( v_1, [ N_x1 , N_x2 ] );
            % u_1 = real(ifft2( v12 )); 
            % v_n(:,2) = v_1;
            % u_n(:,2) = u_1;
            for  i = 2:Ntime

                D1vec = 1i*2*pi*k1vecN;
                D2vec = 1i*2*pi*k2vecN;

                K = k1_l.^2 + ( L_x2/L_x1 * k2_l.^2);
                kvecl2 = K(:);
                
                c = 1 + 1/L_x1;
                
                xi = 3/2 + c*dt - dt*(kvecl2) + L_x1*(kvecl2).^2 ;

                v_02 = reshape( v_0, [ N_x1 , N_x2 ] );
                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                w0 = multiply2D( v_02 , v_02 , 'fourier' );
                w1 = multiply2D( v_12 , v_12 , 'fourier' );
                A0 = (-1/2) * ( D1vec .* w0(:) );
                B0 = (-1/2) * ( D2vec .* w0(:) );
                A1 = (-1/2) * ( D1vec .* w1(:) );
                B1 = (-1/2) * ( D2vec .* w1(:) );

                v_2 = (1./xi) .* ( ...
                    ( ( 2 + ( 2 * c * dt ) ) * v_1 ) - ( 1/2 + ( c * dt ) * v_0 ) + ...
                    ( 2 * dt * ( A1 +  (L_x2/L_x1 * B1) ) ) - ...
                    ( dt * ( A0 +  (L_x2/L_x1 * B0) ) ) ...
                    );

                v_02 = reshape( v_1, [ N_x1 , N_x2 ] );
                u_1 = real(ifft2(v_02));  
                u_1(end,:) = u_1(1,:); % boundary conditions
                u_1(:,end) = u_1(:,1); % boundary conditions
                v_02 = fft2(u_1);
                v_0 = v_02(:);

                v_12 = reshape( v_2, [ N_x1 , N_x2 ] );
                u_2 = real(ifft2(v_12));  
                u_2(end,:) = u_2(1,:); % boundary conditions
                u_2(:,end) = u_2(:,1); % boundary conditions
                v_12 = fft2(u_2);
                v_1 = v_12(:);

                if save_each == 1
                    v_n(:,i) = v_2(:);
                    u_n(:,i) = u_2(:);
                elseif mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_2(:);
                    u_n(:,i/save_each + 1) = u_2(:);
                end
                    
            end
        case 'cnrkw3' % Solve equation by CNRKW3 method
            v_0step = v_n(:,1);
            for i = 2:Ntime
                
                % coefficient terms (Le and Moin 1991)
                alpha_I = [ 4/15 , 1/15 , 1/6 ];
                beta_I = [ 4/15 , 1/15 , 1/6 ];
                alpha_E = [ 8/15 , 5/12 , 3/4 ];
                beta_E = [ 0 , -17/60 , -5/12 ];

                % fourier derivative multiplier
                D1vec = 1i*2*pi*k1vecN;
                D2vec = 1i*2*pi*k2vecN;

                % linear terms
                % Lin = (-1) * ( I^2 * (k1_1.^2 + k2_1.^2) ... 
                    % + I^4 * (k1_1.^4 + 2*(k1_1.^4 .* k2_1.^2) + k2_1.^4));  
                % linear terms
                Lin = (-1) * ( 1i^2 * (kvecl2) + 1i^4 * (kvecl4) );  

                %{
                % nonlinear terms and solution substeps
                Nonlin_v0 = (-1/2) * ...
                    ( D1 .* multiply2D( v_0step , v_0step , 'fourier' ) + ...
                    D2 .* multiply2D( v_0step , v_0step , 'fourier' ) );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );

                Nonlin_v1 = (-1/2) * ...
                    ( D1 .* multiply2D( v_s1 , v_s1 , 'fourier' ) + ...
                    D2 .* multiply2D( v_s1 , v_s1 , 'fourier' ) );
                v_s2 = ( 1 - (dt * alpha_I(2) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(2) * Lin) ) .* v_s1 + ... 
                    (dt * alpha_E(2) * Nonlin_v1) + ...
                    (dt * beta_E(2) * Nonlin_v0) );

                Nonlin_v2 = (-1/2) * ...
                    ( D1 .* multiply2D( v_s2 , v_s2 , 'fourier' ) + ...
                    D2 .* multiply2D( v_s2 , v_s2 , 'fourier' ) );
                v_1 = ( 1 - (dt * alpha_I(3) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(3) * Lin) ) .* v_s2 + ...
                    (dt * alpha_E(3) * Nonlin_v2) + ...
                    (dt * beta_E(3) * Nonlin_v1) ); 
                %}

                %
                % nonlinear terms and solution substeps
                v_0step2 = reshape( v_0step, [ N_x1 , N_x2 ] );
                w = multiply2D( v_0step2 , v_0step2 , 'fourier' );
                Nonlin_v0 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );

                v_s12 = reshape( v_s1, [ N_x1 , N_x2 ] );
                w = multiply2D( v_s12 , v_s12 , 'fourier' );
                Nonlin_v1 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_s2 = ( 1 - (dt * alpha_I(2) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(2) * Lin) ) .* v_s1 + ... 
                    (dt * alpha_E(2) * Nonlin_v1) + ...
                    (dt * beta_E(2) * Nonlin_v0) );

                v_s22 = reshape( v_s2, [ N_x1 , N_x2 ] );
                w = multiply2D( v_s22 , v_s22 , 'fourier' );
                Nonlin_v2 = (-1/2) * ( D1vec .* w(:) + D2vec .* w(:) );
                v_1 = ( 1 - (dt * alpha_I(3) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(3) * Lin) ) .* v_s2 + ...
                    (dt * alpha_E(3) * Nonlin_v2) + ...
                    (dt * beta_E(3) * Nonlin_v1) ); 
                %}

                v_0step = v_1;
                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                u_1 = real(ifft2(v_12));  

                if mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = v_12(:);
                    u_n(:,i/save_each + 1) = ((-1)^(2))*u_1(:);
                end

            end
        case 'imexrk4' % Solve equation by IMEXRK4 method
            v_0step = v_0;
            for i = 2:Ntime
                
                % coefficient terms (Alimo et al. 2021, Table 2)
                alpha_I = [ 343038331393/1130875731271 , 288176579239/1140253497719 ,...
                    253330171251/677500478386 , 189462239225/1091147436423 ];
                beta_I = [ 35965327958/140127563663 , 19632212512/2700543775099 ,...
                    -173747147147/351772688865 , 91958533623/727726057489 ];
                alpha_E = [ 14/25 , 777974228744/1346157007247 ,...
                    251277807242/1103637129625 , 113091689455/220187950967 ];
                beta_E = [ 0 , -251352885992/790610919619 ,...
                    -383714262797/1103637129625 , -403360439203/1888264787188 ];

                % fourier derivative multiplier
                I = 1i;

                % linear terms
                Lin = (-1) * ( I^2 * (k1_l.^2 + k2_l.^2) ... 
                    + I^4 * (k1_l.^4 + k2_l.^4));  

                % nonlinear terms and solution substeps
                Nonlin_v0 = (-1) * (I/2) * ...
                    ( k1_n .* multiply2D(v_0step,v_0step,'fourier') + ...
                    k2_n .* multiply2D(v_0step,v_0step,'fourier') );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );

                Nonlin_v1 = (-1) * (I/2) * ...
                    ( k1_n .* multiply2D(v_s1,v_s1,'fourier') + ...
                    k2_n .* multiply2D(v_s1,v_s1,'fourier') );
                v_s2 = ( 1 - (dt * alpha_I(2) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(2) * Lin) ) .* v_s1 + ... 
                    (dt * alpha_E(2) * Nonlin_v1) + ...
                    (dt * beta_E(2) * Nonlin_v0) );

                Nonlin_v2 = (-1) * (I/2) * ...
                    ( k1_n .* multiply2D(v_s2,v_s2,'fourier') + ...
                    k2_n .* multiply2D(v_s2,v_s2,'fourier') );
                v_s3 = ( 1 - (dt * alpha_I(3) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(3) * Lin) ) .* v_s2 + ...
                    (dt * alpha_E(3) * Nonlin_v2) + ...
                    (dt * beta_E(3) * Nonlin_v1) ); 

                Nonlin_v3 = (-1) * (I/2) * ...
                    ( k1_n .* multiply2D(v_s3,v_s3,'fourier') + ...
                    k2_n .* multiply2D(v_s3,v_s3,'fourier') );
                v_1 = ( 1 - (dt * alpha_I(4) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(4) * Lin) ) .* v_s3 + ...
                    (dt * alpha_E(4) * Nonlin_v3) + ...
                    (dt * beta_E(4) * Nonlin_v2) ); 

                v_0step = v_1;
                u_1 = real(ifft2(v_1));  

                if mod(i,save_each) == 0
                    v_n(:,i/save_each + 1) = reshape( v_1 , [ N_x1 * N_x2 , 1 ] );
                    u_n(:,i/save_each + 1) = reshape( ((-1)^(i-1))*u_1 , [ N_x1 * N_x2 , 1 ] );
                end

            end
        
        case 'etdrk4' % Solve equation by mETDRK4 method
            %%%% under construction
            %{
            % v_0step = v_0;
            for i=2:Ntime
                
                % fourier derivative multiplier
                I = 2*pi*1i;

                % linear terms
                Lin = (-1) * ( I^2 * (k1_1.^2 + k2_1.^2) ... 
                    + I^4 * (k1_1.^4 + k2_1.^4));  

                % Precompute various ETDRK4 scalar quantities:
                h = dt; % time step
                L = Lin; % Fourier multipliers
                E = exp(h*L); 
                E2 = exp(h*L/2);
                M = 16; % no. of points for complex means
                r = exp(1i*pi*((1:M)-.5)/M); % roots of unity
                LR = h*L(:,ones(M,1)) + r(ones(N_x1,1),:);
                Q = h*real(mean( (exp(LR/2)-1)./LR ,2));
                f1 = h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
                f2 = h*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
                f3 = h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
                
                %{
                % nonlinear terms and solution substeps
                Nonlin_v0 = (-1) * (I/2) * ...
                    ( k1_0 .* multiply2D(v_0step,v_0step,'fourier') + ...
                    k2_0 .* multiply2D(v_0step,v_0step,'fourier') );
                v_s1 = ( 1 - (dt * alpha_I(1) * Lin) ).^(-1) .* ...
                    ( ( 1 + (dt * beta_I(1) * Lin) ) .* v_0step + ...
                    (dt * alpha_E(1) * Nonlin_v0) + ...
                    (dt * beta_E(1)) );
                %}

                % Main time-stepping loop:
                uu = u_0(:,1); 
                v = v_0(:,1);
                tt = 0;
                tmax = 150; 
                nmax = round(tmax/h); 
                nplt = floor((tmax/100)/h);
                g = -0.5i*k1_0(1,:)';
                for n = 1:nmax
                    t = n*h;
                    Nv = g.*fft(real(ifft(v)).^2);
                    a = E2.*v + Q.*Nv;
                    Na = g.*fft(real(ifft(a)).^2);
                    b = E2.*v + Q.*Na;
                    Nb = g.*fft(real(ifft(b)).^2);
                    c = E2.*a + Q.*(2*Nb-Nv);
                    Nc = g.*fft(real(ifft(c)).^2);
                    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
                    if mod(n,nplt)==0
                        u = real(ifft(v));
                        uu = [uu,u]; 
                        tt = [tt,t];
                    end
                end

                %{ 
                % use this idea for innermost loop above (uu, tt)
                x1s = size( x1 , size(x1,2) * strips );
                x2s = size( x2 , size(x2,2) * strips );
                x1s( x1 , 1:size(x1,2) ) = x1;
                x2s( x2 , 1:size(x2,2) ) = x2;
                for i = 2:strips
                    col_start_x1 = (i-1) * size(x1,2) + 1;
                    col_end_x1 = (i) * size(x1,2);
                    col_start_x2 = (i-1) * size(x2,2) + 1;
                    col_end_x2 = (i) * size(x2,2);
                    x1s( x1 , col_start_x1:col_end_x1 ) = [ x1s ; x1 ];
                    x2s( x2 , col_start_x2:col_end_x2 ) = [ x2s ; x2 ];
                end
                %}

                % Plot results:
                surf(tt,x1(1,:),uu), shading interp, lighting phong, axis tight
                view([-90 90]), colormap(autumn); 
                set(gca,'zlim',[-5 50])
                light('color',[1 1 0],'position',[-1,2,2])
                material([0.30 0.60 0.60 40.00 1.00]);
                                    
            end
            %}
    end
end
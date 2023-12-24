function [ v_n , u_n ] = DirectSolve_2DKS_v2_0(IC,method,N_x1,L_s1,L_s2,dt,T,save_each)

%{ 
Output:
v_n: solution vector for each time step in Fourier space
u_n: solution vector for each time step in Physical space
%}

%%% (1) initialize problem %%%

    % Discretized equally in each dimension
    N_x2 = N_x1;

    % length-scale parameters
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
    K = k1_l.^2 + k2_l.^2; % 2nd order linear term
    kvecl2 = K(:);
    KK = k1_l.^4 + k2_l.^4;
    kvecl4 = KK(:);
    
    switch method
        case 'imexrk4vec2a'
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
        Ntime_save = 1;
        Ntime = time1;
        save_each = time1;
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
    % figure(1); 
    % surf(k1_0,k2_0,abs(v_i)); set(gca,'xscale','log','yscale','log','zscale','log');
    % view(90,0);
    figure(2);
    surfc(x1,x2,u_i);
    shading interp
    pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );
    view(3);
    %}

    switch method 

        case 'imexrk4vec2a' % Solve vectorized equation by IMEXRK4 method
            v_0step = v_n(:,1);
            for i = 2:Ntime

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

                v_12 = reshape( v_1, [ N_x1 , N_x2 ] );
                u_1 = real(ifft2(v_12));  
                v_0step = v_1;

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
function [maxL2inT,u_IC,u_TC,energy,v_mean,projcoeffradialevolution,projcoeffmodeevolution] = ...
    process_energy(u_IC, savedata, IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
    parameterlist,optparameters,parfiglist,optparfiglist)

    %% Note: cheaper to recompute rather than read from file!

    %% setup
    
    % time window
    timewindow = linspace(0,T,Ntime);
    
    % length-scale parameters
    L1 = 2*pi*L_s1;
    L2 = 2*pi*L_s2;
    
    % unit physical space domain
    x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
    x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
    [ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
    %% solution grid
    % fourier space domain for nonlinear term
    k1_n_pts = 2*pi/L1*[ 0 : N/2-1 , 0 , -N/2+1 : -1]; 
    k2_n_pts = 2*pi/L2*[ 0 : N/2-1 , 0 , -N/2+1 : -1]; 
    [ k1_n , k2_n ] = meshgrid(k1_n_pts,k2_n_pts);              % 2-dimensional grid

    % fourier space domain for linear term
    k1_l_pts = 2*pi/L1*[0:N/2 -N/2+1:-1];
    k2_l_pts = 2*pi/L2*[0:N/2 -N/2+1:-1];
    [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts);              % 2-dimensional grid
    
    % Fourier space domain
    kx = (2*pi/L1) * (-N/2 : N/2-1);   
    ky = (2*pi/L2) * (-N/2 : N/2-1);
    [KX, KY] = meshgrid(kx, ky);       
    K2 = KX.^2 + KY.^2;                % |k|^2

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
    Bilap = 1i^4 * (kvecl4);                                    % Bi-Laplacian operator
    D1vec = 1i*k1vecN;                                          % Differential operator
    D2vec = 1i*k2vecN;                                          % Differential operator
    
    % coefficient terms (Alimo et al. 2021, Table 2)
    alpha_I = [ 343038331393/1130875731271 , 288176579239/1140253497719 , 253330171251/677500478386 , 189462239225/1091147436423 ];
    beta_I = [ 35965327958/140127563663 , 19632212512/2700543775099 , -173747147147/351772688865 , 91958533623/727726057489 ];
    alpha_E = [ 14/25 , 777974228744/1346157007247 , 251277807242/1103637129625 , 113091689455/220187950967 ];
    beta_E = [ 0 , -251352885992/790610919619 , -383714262797/1103637129625 , -403360439203/1888264787188 ];
    
    % impose initial and boundary conditions in physical and fourier space
    u_0 = reshape( u_IC, [ N , N ]);
    v_0 = fft2(u_0);            % FFT of physical initial condition
    v_step = v_0(:);

    % number of timesteps
    save_each = 1;
    time1 = ceil(T/dt); 
    time2 = ceil(time1/save_each);
    Ntime = max(time1,time2);
    Ntime_save = min(time1,time2);
    save_each = ceil(Ntime/Ntime_save);
    Ntime_save = ceil(Ntime/save_each);
    Ntime = save_each*Ntime_save;
    %%    

    %% process

    % compute L2 energy and Fourier mode evolution
    sol_samples3 = NaN(Ntime,8*3);
    mean_sum_s1 = 0;
    mean_sum_s2 = 0;
    mean_sum_s3 = 0;
    energyL2 = NaN(Ntime,1);
    energyH1 = energyL2;
    energyH2 = energyL2;
    energyL2_lap = energyL2;
    energyL2_bilap = energyL2;
    astripwidth = NaN(Ntime,2);                               % analyticity strip width
    v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1, Ntime+1 ); % uniform-distance radial energy spectrum
    v_mean(:,1) = NaN;
    projcoeffradialevolution = v_mean;
    projcoeffmodeevolution = [];
    Ntime_remaining = Ntime;
    for i = 1:Ntime

        % nonlinear terms and solution substeps
        if i > 1
            [v_1,Nonlin_v1,v_1_t] = ks_nonlinear_fwd( v_step, D1vec, D2vec, N, Lin, dt, alpha_I, beta_I, alpha_E, beta_E );
        else
            v_1 = v_step;
            v_02 = reshape( v_1, [ N , N ] );   % old term
            %u_nm1 = real(ifft2(v_02));          % old term
        end

        % correction of non-zero mean solution
        if abs(mean(v_1)) > 1e-5 
            v_1(1) = 0;
        end

        %%% start equation check %%%

        % save solution step to workspace
        v_12 = reshape( v_1, [ N , N ] );   % new term
        u_n = real(ifft2(v_12));            % new term
        v_step = v_1;                       % new term

        % compute 2D KS terms
        if i > 1
            %v_step2x 	= reshape( D1vec .* v_step, [ N , N ] ); % f_x 
            %v_step2y 	= reshape( D2vec .* v_step, [ N , N ] ); % f_y 
            v_2lap 		= reshape( Lap .* v_step, [ N , N ] ); 	 % lap(f) 
            v_2bilap 	= reshape( Bilap .* v_step, [ N , N ] ); % bilap(f) 
            %{
            w1_r 		= multiply2D( v_step2x , v_step2x , 'fourier2real' ); 	% f_x * f_x in physical space (pseudospectral) 
            w1s_r 		= multiply2D( v_step2y , v_step2y , 'fourier2real' ); 	% f_y * f_y in physical space (pseudospectral) 
            Nonlin_v1_r = (1/2) * ( w1_r + w1s_r ); 				            % (1/2)*(f_x * f_x + f_y * f_y) in physical space 
            Nonlin_v1 	= multiply2D( fft2(Nonlin_v1_r) , 0 ,'dealias'); 	    % 2/3 dealias 
            %}
    
            %u_t         = (u_n - u_nm1) / dt;           % physical time derivative
            u_t 	    = real(ifft2(v_1_t)); 		    % physical time derivative
            u_nonlin 	= real(ifft2(Nonlin_v1)); 		% physical nonlinear term
            u_lap 		= real(ifft2(v_2lap)); 			% physical laplacian term
            u_bilap 	= real(ifft2(v_2bilap));		% physical bilaplacian term
            sum_u       = u_t + u_nonlin + u_lap + u_bilap; 
            sumf = fft2(sum_u); 
            sumf(1,1) = 0; 
            sum_u = real(ifft2(sumf));                  % check mean-zero sum equal to 0
    
            % sample indices
            if i == 2
                [~,maxidx] = max(abs(u_n(:)));
                [i1, j1] = ind2sub(size(u_n), maxidx);
                i2 = 1; j2 = 1;
                i3 = ceil(N/2); j3 = ceil(N/2);
            end
    
            u_n_s1      = u_n(i1,j1);       % solution sample 1
            u_t_s1      = u_t(i1,j1);       % physical time derivative sample 1
            u_nonlin_s1 = u_nonlin(i1,j1);  % physical nonlinear term sample 1
            u_lap_s1    = u_lap(i1,j1); 	% physical laplacian term sample 1
            u_bilap_s1 	= u_bilap(i1,j1);	% physical bilaplacian term sample 1
            sum_s1      = sum_u(i1,j1);     % check sum equal to 0
            mean_sum_s1 = (mean_sum_s1 + sum_s1)/i; % check average sum equal to 0
            rel_s1 = abs(sum_s1) / (abs(u_t_s1) + abs(u_nonlin_s1) + abs(u_lap_s1) + abs(u_bilap_s1) );
    
            u_n_s2      = u_n(i2,j2);       % solution sample 2
            u_t_s2      = u_t(i2,j2);       % physical time derivative sample 2
            u_nonlin_s2 = u_nonlin(i2,j2);  % physical nonlinear term sample 2
            u_lap_s2    = u_lap(i2,j2); 	% physical laplacian term sample 2
            u_bilap_s2 	= u_bilap(i2,j2);	% physical bilaplacian term sample 2
            sum_s2      = sum_u(i2,j2);     % check sum equal to 0
            mean_sum_s2 = (mean_sum_s2 + sum_s2)/i; % check average sum equal to 0
            rel_s2 = abs(sum_s2) / (abs(u_t_s2) + abs(u_nonlin_s2) + abs(u_lap_s2) + abs(u_bilap_s2) );
    
            u_n_s3      = u_n(i3,j3);       % solution sample 3
            u_t_s3      = u_t(i3,j3);       % physical time derivative sample 3
            u_nonlin_s3 = u_nonlin(i3,j3);  % physical nonlinear term sample 3
            u_lap_s3    = u_lap(i3,j3); 	% physical laplacian term sample 3
            u_bilap_s3 	= u_bilap(i3,j3);	% physical bilaplacian term sample 3
            sum_s3      = sum_u(i3,j3);     % check sum equal to 0
            mean_sum_s3 = (mean_sum_s3 + sum_s3)/i; % check average sum equal to 0
            rel_s3 = abs(sum_s3) / (abs(u_t_s3) + abs(u_nonlin_s3) + abs(u_lap_s3) + abs(u_bilap_s3) );
    
            sol_samples3(i,:) = [ u_n_s1 , u_n_s2 , u_n_s3 , ...
                u_t_s1 , u_t_s2 , u_t_s3 , ...
                u_nonlin_s1 , u_nonlin_s2 , u_nonlin_s3 , ...
                u_lap_s1 , u_lap_s2 , u_lap_s3 , ...
                u_bilap_s1 , u_bilap_s2 , u_bilap_s3 , ...
                sum_s1 , sum_s2 , sum_s3 , ...
                mean_sum_s1 , mean_sum_s2 , mean_sum_s3 , ...
                rel_s1 , rel_s2 , rel_s3  ];
    
            %u_nm1 = u_n;
        end
        %%% end equation check %%%



        % save solution step to workspace
        v_12 = reshape( v_1, [ N , N ] );   % new term
        u_n = real(ifft2(v_12));            % new term
        v_step = v_1; 

        %{
        if Ntime < Ntime_save_max && i == 1
            [u_n, v_n] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime T], utility1);
            %u_n =  real(ifft2(ifftshift(v_n)));
            fprintf('Computed energy evolution at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
        else
            if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                currentT = (i+Ntime_save_max-1)/Ntime*T;
                [u_n, v_n] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, [Ntime_save_max T], utility1);
                %u_n =  real(ifft2(ifftshift(v_n)));
                Ntime_remaining = Ntime_remaining - Ntime_save_max;
            elseif (mod(i,Ntime_save_max) == 1)
                [u_n, v_n] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime_remaining T], utility1);
                %u_n =  real(ifft2(ifftshift(v_n)));
                fprintf('Computed energy evolution at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            end
        end

        if mod(i,Ntime_save_max) ~= 0
            imod = mod(i,Ntime_save_max);
        else
            imod = Ntime_save_max;
        end
        %}

        % projection coefficients and unstable modes
        [~,projcoeffs,unstablemodes] = eigenfunction_validation(u_n(:),L_s1, N, T,IC,'full');
        if size(projcoeffmodeevolution,1) ~= size(unstablemodes,1)
            projcoeffmodeevolution = [projcoeffmodeevolution ; zeros(size(unstablemodes,1)-size(projcoeffmodeevolution,1),Ntime+2)];
        end
        for j = 1:size(unstablemodes,1)
            % radial projection sum
            radius = sqrt((unstablemodes(j,1))^2+(unstablemodes(j,2))^2);
            if ~isempty(find(abs(projcoeffradialevolution(:,1) - radius) < 1e-10 , 1))
                idx = find(abs(projcoeffradialevolution(:,1) - radius) < 1e-10 , 1);
                projcoeffradialevolution(idx,i+1) = projcoeffradialevolution(idx,i+1) + projcoeffs(j).^2/norm(projcoeffs).^2;%projcoeffs(j,1);
            else
                idx = find(isnan(projcoeffradialevolution(:,1)) , 1, 'first');
                projcoeffradialevolution(idx,1) = radius;
                projcoeffradialevolution(idx,i+1) = projcoeffs(j).^2/norm(projcoeffs).^2;%projcoeffs(j,1);
            end

            % mode projection coefficient
            if ~isempty(find(abs(projcoeffmodeevolution(:,1) - unstablemodes(j,1)) < 1e-10 & abs(projcoeffmodeevolution(:,2) - unstablemodes(j,2)) < 1e-10, 1, 'first'))
                idx = find(abs(projcoeffmodeevolution(:,1) - unstablemodes(j,1)) < 1e-10 & abs(projcoeffmodeevolution(:,2) - unstablemodes(j,2)) < 1e-10, 1, 'first');
                projcoeffmodeevolution(idx,i+2) = projcoeffs(j).^2/norm(projcoeffs).^2;
            else
                idx = find(projcoeffmodeevolution(:,1) == 0 & projcoeffmodeevolution(:,2) == 0, 1, 'first');
                projcoeffmodeevolution(idx,1:2) = unstablemodes(j,1:2);
                projcoeffmodeevolution(idx,i+2) = projcoeffs(j).^2/norm(projcoeffs).^2;
            end

        end

        u_i = reshape( u_n(:) , [ N , N ] );
        v = fftshift(abs(fft2(u_i)).^2);
        %energyL2(i,1) = (sum( u_n(:) .* conj(u_n(:)) )*(L1*L2)/(N*N));
        energyL2(i,1) = sum( v(:) )*(L1*L2)/(N*N)^2;
        energyH1(i,1) = sum( (1 + K2(:)) .* v(:) )*(L1*L2)/(N*N)^2;
        energyH2(i,1) = sum( (1 + K2(:)).^2 .* v(:) )*(L1*L2)/(N*N)^2;
        energyL2_lap(i,1) = sum( (K2(:)).^2 .* v(:) )*(L1*L2)/(N*N)^2;
        energyL2_bilap(i,1) = sum( (K2(:)).^4 .* v(:) )*(L1*L2)/(N*N)^2;

        for j = 1:N
            for k = 1:N
                radius = sqrt((j-(N/2+1))^2+(k-(N/2+1))^2);
                if ~isempty(find(abs(v_mean(:,1) - round(radius)) < 1e-10 , 1))
                    idx = find(abs(v_mean(:,1) - round(radius)) < 1e-10 , 1);
                    v_mean(idx,i+1) = v_mean(idx,i+1) + v(j,k);
                else
                    idx = find(isnan(v_mean(:,1)) , 1, 'first');
                    v_mean(idx,1) = round(radius);
                    v_mean(idx,i+1) = v(j,k);
                end
            end
        end
        v_mean = sortrows(v_mean,1);
        v_mean(:,i+1) = v_mean(:,i+1)./2;
            
        % Width of the analyticity strip [C , delta]
        if i > 1
            [astripwidth(i,1),astripwidth(i,2)] = expfit_delta(v_mean(2:end,1), v_mean(2:end,i+1), astripwidth(i-1,2));
        else
            [astripwidth(i,1),astripwidth(i,2)] = expfit_delta(v_mean(2:end,1), v_mean(2:end,i+1), 1e15);
        end

        if i == 1
            u_IC = u_n;
        elseif i == Ntime
            u_TC = u_n;
        end

    end
    v_mean = v_mean(2:end,:);
    
    %{
    % L2 energy time derivative computation
    energyL2_t = NaN(Ntime-1,1);
    dt_save = T/(Ntime-1);
    for i = 2:length(energyL2)
        energyL2_t(i-1,1) = ( energyL2(i,1) - energyL2(i-1,1) ) / dt_save;
    end
    %}

    % max energy in time window
    [ maxL2, maxL2index ] = max(energyL2);
    maxL2time = timewindow(maxL2index);
    maxL2inT = [ maxL2time , maxL2 ];

    energy = [ timewindow', energyL2 , energyH1 , energyH2, energyL2_lap, ...
        astripwidth(:,1) , astripwidth(:,2), energyL2_bilap, sol_samples3 ];

    if savedata == 1
        energydata_file = [pwd '/data/energy/energy_' parameterlist '.dat'];
        spectrum_file = [pwd '/data/spectrum/spectrum_' parameterlist '.dat'];
        switch IC
            case {'optimized'}
                energydata_file = [pwd '/data/energy/energy_' optparameters '.dat'];
                spectrum_file = [pwd '/data/spectrum/spectrum_' optparameters '.dat'];
        end
        writematrix(energy, energydata_file,'Delimiter','tab');
        writematrix(v_mean, spectrum_file,'Delimiter','tab');
    end
end
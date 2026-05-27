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
    L_ref = 1.02;
    fraclist = [];%[ 1/2, 1/3, 1/4, 1/6, 1/12, 1/16 ];
    samplecount = length(fraclist) + 1;
    sol_samplemax = NaN(Ntime,8*(samplecount+1));
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
            v_step2x 	= reshape( D1vec .* v_step, [ N , N ] );    % f_x 
            v_step2y 	= reshape( D2vec .* v_step, [ N , N ] );    % f_y 
            v_2lap 		= reshape( Lap .* v_step, [ N , N ] ); 	    % lap(f) 
            v_2bilap 	= reshape( Bilap .* v_step, [ N , N ] );    % bilap(f) 
            u_t 	    = real(ifft2(v_1_t)); 		                % physical time derivative
            u_nonlin 	= real(ifft2(Nonlin_v1)); 		            % physical nonlinear term
            u_lap 		= real(ifft2(v_2lap)); 			            % physical laplacian term
            u_bilap 	= real(ifft2(v_2bilap));		            % physical bilaplacian term
            u_x         = real(ifft2(v_step2x));		            % physical dim 1 derivative
            u_y         = real(ifft2(v_step2y));		            % physical dim 2 derivative
            sum_u       = u_t + u_nonlin + u_lap + u_bilap;         % residual of solution
            sumf        = fft2(sum_u); 
            sumf(1,1)   = 0;                                        % mean-zero correction of residual
            sum_u       = real(ifft2(sumf));                        % check mean-zero residual equal to 0
    
            % sample indices
            if i == 2

                idxlist = NaN(samplecount,2);

                [~,maxidx] = max(abs(u_n(:)));
                [idxlist(1,1), idxlist(1,2)] = ind2sub(size(u_n), maxidx);

                for ii = 2:samplecount
                    fracidx = fraclist(ii-1);
                    geoshift = round(N*fracidx*(L_ref/L_s1)*(L_ref/L_s2));
                    if startsWith(optparameters,'tg')
                        idxlist(ii,1) = max( mod(idxlist(1,1) + geoshift,N) , 1);
                    else
                        idxlist(ii,1) = idxlist(1,1);
                    end
                    idxlist(ii,2) = max( mod(idxlist(1,2) + geoshift,N) , 1);
                end

                row = idxlist(1,1);
                col = idxlist(1,2);
                sz = size(u_n);
                
                % periodic mean for rows
                theta_r = 2*pi*(row-1)/sz(1);
                rowmin = mod(atan2(mean(sin(theta_r)), mean(cos(theta_r))) ...
                             * sz(1)/(2*pi), sz(1)) + 1;
                
                % periodic mean for columns
                theta_c = 2*pi*(col-1)/sz(2);
                colmin = mod(atan2(mean(sin(theta_c)), mean(cos(theta_c))) ...
                             * sz(2)/(2*pi), sz(2)) + 1;
                
                rowmin = round(rowmin);
                colmin = round(colmin);
                [Ny, Nx] = size(u_n);
                rowc = floor((Ny+1)/2) + 1;
                colc = floor((Nx+1)/2) + 1;
                %drow = rowc - rowmin;
                %dcol = colc - colmin;

            end
    
            for idx = 1:samplecount
                u_n_s1      = u_n(idxlist(idx,1),idxlist(idx,2));       % solution sample 1
                u_t_s1      = u_t(idxlist(idx,1),idxlist(idx,2));       % physical time derivative sample 1
                u_nonlin_s1 = u_nonlin(idxlist(idx,1),idxlist(idx,2));  % physical nonlinear term sample 1
                u_lap_s1    = u_lap(idxlist(idx,1),idxlist(idx,2)); 	% physical laplacian term sample 1
                u_bilap_s1 	= u_bilap(idxlist(idx,1),idxlist(idx,2));	% physical bilaplacian term sample 1
                sum_s1      = sum_u(idxlist(idx,1),idxlist(idx,2));     % check sum equal to 0
                u_x_s1      = u_x(idxlist(idx,1),idxlist(idx,2));       % physical dim 1 derivative sample 1
                u_y_s1      = u_y(idxlist(idx,1),idxlist(idx,2));       % physical dim 2 derivative sample 1
    
                sol_samplemax(i,(idx)*8 + 1:(idx)*8 + 8) = [ u_n_s1 , u_t_s1 , u_nonlin_s1 , u_lap_s1 , u_bilap_s1 , sum_s1 , u_x_s1 , u_y_s1 ];
            end

            %{
            geoshift_ref = round(N*1/2*(L_ref/L_s1)*(L_ref/L_s2));
            Nstart = max(round(N/2 - geoshift_ref + 1),1);
            Nend = min(round(N/2 + geoshift_ref),N);
            if startsWith(optparameters,'tg')
                u_n_ps = circshift(u_n, [drow, dcol]);
                subN = Nend-Nstart+1;
                u_n_ref = zeros(subN);
                x0 = N/2;
                y0 = 0;
                xs = 1;
                for xshift = 1:N/2
                    ys = 1;
                    for yshift = 1:N/2
                        u_n_ref(xs,ys) = u_n_ps(x0+yshift,y0+yshift);
                        ys = ys + 1;
                    end
                    xs = xs + 1;
                    x0 = N/2 - xshift;
                    y0 = 0 + xshift;
                end


            else
                u_n_ps = circshift(u_n, [drow, dcol]);
                u_n_ref = u_n_ps(Nstart:Nend,:);


            end
            %}

            %v_12_ref = fft2(u_n_ref);
            v = fftshift(abs(v_12).^2);
            u_n_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(v_1_t).^2);
            u_t_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(Nonlin_v1).^2);
            u_nonlin_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(v_2lap).^2);
            u_lap_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(v_2bilap).^2);
            u_bilap_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(sumf).^2);
            sum_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(v_step2x).^2);
            u_x_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            v = fftshift(abs(v_step2y).^2);
            u_y_L2 = sum( v(:) )*(L1*L2)/(N*N)^2;

            sol_samplemax(i, 1:8) = [ u_n_L2 , u_t_L2 , u_nonlin_L2 , u_lap_L2 , u_bilap_L2 , sum_L2 , u_x_L2 , u_y_L2 ];

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
        %energyL2(i,1) = (sum( u_n(:) .* conj(u_n(:)) )*(L1*L2)/(N*N));
        v = fftshift(abs(fft2(u_i)).^2);
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
        astripwidth(:,1) , astripwidth(:,2), energyL2_bilap, sol_samplemax ];

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

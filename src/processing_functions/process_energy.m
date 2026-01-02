function [maxL2inT,u_IC,u_TC,energyL2,energyH1,energyH2,astripwidth,v_mean,projcoeffradialevolution,projcoeffmodeevolution] = ...
    process_energy(savedata, IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
    parameterlist,optparameters,parfiglist,optparfiglist)

    %% setup
    
    % time window
    timewindow = linspace(0,T,Ntime);
    
    % length-scale parameters
    L1 = 2*pi*L_s1;
    L2 = 2*pi*L_s2;
    
    % unit physical space domain
    x1_pts = L_s1*linspace( 0 , 1 - 1/N , N ); 
    x2_pts = L_s2*linspace( 0 , 1 - 1/N , N ); 
    [ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
    % unit physical space domain
    x1_pts = 2*pi*L_s1*linspace( 0 , 1 - 1/N , N ); 
    x2_pts = 2*pi*L_s2*linspace( 0 , 1 - 1/N , N ); 
    [ x1pi , x2pi ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
    % Fourier space domain
    kx = (2*pi/L1) * (-N/2 : N/2-1);   
    ky = (2*pi/L2) * (-N/2 : N/2-1);
    [KX, KY] = meshgrid(kx, ky);       
    K2 = KX.^2 + KY.^2;                % |k|^2
    
    % 4-period physical space domain
    x12_pts = 2*L_s1*linspace( 0 , 1 - 1/N , 2*N ); 
    x22_pts = 2*L_s2*linspace( 0 , 1 - 1/N , 2*N ); 
    [ x12x , x22x ] = meshgrid(x12_pts,x22_pts); % 2-dimensional grid

    %% process

    % compute L2 energy and Fourier mode evolution
    energyL2 = NaN(Ntime,1);
    energyH1 = energyL2;
    energyH2 = energyL2;
    astripwidth = NaN(Ntime,2);                               % analyticity strip width
    v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1, Ntime+1 ); % uniform-distance radial energy spectrum
    v_mean(:,1) = NaN;
    projcoeffradialevolution = v_mean;
    projcoeffmodeevolution = [];
    Ntime_remaining = Ntime;
    for i = 1:Ntime

        if Ntime < Ntime_save_max && i == 1
            [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime T], utility1);
            fprintf('Computed energy evolution at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
        else
            if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                currentT = (i+Ntime_save_max-1)/Ntime*T;
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, [Ntime_save_max T], utility1);
                Ntime_remaining = Ntime_remaining - Ntime_save_max;
            elseif (mod(i,Ntime_save_max) == 1)
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime_remaining T], utility1);
                fprintf('Computed energy evolution at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))
            end
        end
        
        if mod(i,Ntime_save_max) ~= 0
            imod = mod(i,Ntime_save_max);
        else
            imod = Ntime_save_max;
        end

        % projection coefficients and unstable modes
        [~,projcoeffs,unstablemodes] = eigenfunction_validation(u_n(:,imod),L_s1, N, T,IC,'full');
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

        u_i = reshape( u_n(:,imod) , [ N , N ] );
        v = fftshift(abs(fft2(u_i)).^2);
        %energyL2(i,1) = (sum( u_n(:,imod) .* conj(u_n(:,imod)) )*(L1*L2)/(N*N));
        energyL2(i,1) = sum( v(:) )*(L1*L2)/(N*N)^2;
        energyH1(i,1) = sum( (1 + K2(:)) .* v(:) )*(L1*L2)/(N*N)^2;
        energyH2(i,1) = sum( (1 + K2(:)).^2 .* v(:) )*(L1*L2)/(N*N)^2;

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
            u_IC = u_n(:,1);
        elseif i == Ntime
            u_TC = u_n(:,end);
        end

    end
    v_mean = v_mean(2:end,:);
    
    % L2 energy time derivative computation
    energyL2_t = NaN(Ntime-1,1);
    dt_save = T/(Ntime-1);
    for i = 2:length(energyL2)
        energyL2_t(i-1,1) = ( energyL2(i,1) - energyL2(i-1,1) ) / dt_save;
    end

    % max energy in time window
    [ maxL2, maxL2index ] = max(energyL2);
    maxL2time = timewindow(maxL2index);
    maxL2inT = [ maxL2time , maxL2 ];
    
    if savedata == 1
        energyL2data_file = [pwd '/data/energyL2/energyL2_' parameterlist '.dat'];
        energyL2deriv_file = [pwd '/data/energyL2_t/energyL2_t_' parameterlist '.dat'];
        spectrum_file = [pwd '/data/spectrum/spectrum_' parameterlist '.dat'];
        switch IC
            case {'optimized'}
                energyL2data_file = [pwd '/data/energyL2/energyL2_' optparameters '.dat'];
                energyL2deriv_file = [pwd '/data/energyL2_t/energyL2_t_' optparameters '.dat'];
                spectrum_file = [pwd '/data/spectrum/spectrum_' optparameters '.dat'];
        end
        writematrix(energyL2, energyL2data_file,'Delimiter','tab');
        writematrix(energyL2_t, energyL2deriv_file,'Delimiter','tab');
        writematrix(v_mean, spectrum_file,'Delimiter','tab');
    end
end
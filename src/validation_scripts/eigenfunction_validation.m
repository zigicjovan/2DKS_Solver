function [match_score,ampstars,modes] = eigenfunction_validation(u_IC_opt, L_s1, N, T, IC, type)

    L_s2 = L_s1;
    
    u = reshape(u_IC_opt,[N,N]);
    %[~, max_idx] = max(u(:));
    [~, max_idx] = min(u(:));
    [rowmax, colmax] = ind2sub(size(u), max_idx);
    u1 = [ u(rowmax+1:end,:) ;  u(1:rowmax,:)  ];
    u2 = [ u1(:,colmax:end) ,  u1(:,1:colmax-1)  ];
    
    %u2 = u;

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
    
    x1_pts = L_s1*linspace( 0 , 1 - 1/N , N ); 
    x2_pts = L_s2*linspace( 0 , 1 - 1/N , N ); 
    [ x1p , x2p ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
    amps = ceil(L_s1)^2;
    modes = NaN(amps,2);
    modelist = NaN(amps,2);
    k1c = 0;
    k2c = 1;
    for k = 1:size(modelist,1)
        if isnan(modelist(k,1))
            if k1c == 0
                modelist(k,:) = [ k1c k2c ];
                modelist(k+1,:) = flip(modelist(k,:));
                k1c = k1c + 1;
            elseif k2c > k1c
                modelist(k,:) = [ k1c k2c ];
                modelist(k+1,:) = flip(modelist(k,:));
                k1c = k1c + 1;
            elseif k2c == k1c
                modelist(k,:) = [ k1c k2c ];
                k2c = k2c + 1;
                k1c = 0;
            end
        end
    end

    modelist = [ modelist ; -modelist];

    for k = 1:length(modelist)/2
        if modelist(k,1) > 0 && modelist(k,2) > 0
            modelist = [ modelist ; -modelist(k,1),modelist(k,2) ; modelist(k,1),-modelist(k,2) ];
        end
    end

    modecount = 1;
    switch type
        case {'active','full'}
            for k = 1:size(modelist,1)
                if sqrt(modelist(k,1)^2 + modelist(k,2)^2) < L_s1
                    modes(modecount,:) = modelist(k,:);
                    modecount = modecount + 1;
                end
            end
        case 'dominant'
            for k = 1:size(modelist,1)
                if abs(L_s1 - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 
                    modes(modecount,:) = modelist(k,:);
                    modecount = modecount + 1;
                end
            end
            if abs(L_s1 - sqrt(3)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(2) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(4) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(6)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(4) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(8) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(9)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(8) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(10) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(13)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(10) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(16) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(17)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(16) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(18) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(19)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(18) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(20) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(23)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(20) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(26) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
            if abs(L_s1 - sqrt(29)) < 1e-10
                for k = 1:size(modelist,1)
                    if abs(sqrt(26) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2 || abs(sqrt(32) - sqrt(modelist(k,1)^2 + modelist(k,2)^2)*sqrt(2)) < 1e-2
                        modes(modecount,:) = modelist(k,:);
                        modecount = modecount + 1;
                    end
                end
            end
    end
    modes = rmmissing(modes);
    amps = size(modes,1);
    ampstars = NaN(amps,1);

    for k = 1:size(modes,1)
        u_ef = cos( (modes(k,1)*L_x1*x1 + modes(k,2)*L_x2*x2) );
        u_ef = u_ef / sqrt(sum((u_ef(:)) .* conj((u_ef(:))) )*(L1*L2)/N^2);
        %[~, max_idx] = max(u_ef(:));
        [~, max_idx] = min(u_ef(:));
        [rowmax, colmax] = ind2sub(size(u_ef), max_idx);
        u1_ef = [ u_ef(rowmax+1:end,:) ;  u_ef(1:rowmax,:)  ];
        u2_ef = [ u1_ef(:,colmax:end) ,  u1_ef(:,1:colmax-1)  ];

        u_centered = u2 - mean(u2(:));
        u_norm = u_centered / sqrt( sum((u_centered(:)) .* conj((u_centered(:))))*(L1*L2)/N^2 );
        phi_norm = u2_ef / sqrt( sum((u2_ef(:)) .* conj((u2_ef(:))))*(L1*L2)/N^2 );
        %phi_norm = u_ef / sqrt( sum((u_ef(:)) .* conj((u_ef(:))))*(L1*L2)/N^2 );
        ampstars(k,1) = (sum((u_norm(:)) .* conj((phi_norm(:))) )*(L1*L2)/N^2);
    end

    % check orthogonality
    orthocheck = NaN(amps);
    for m = 1:size(modes,1)
        mm = cos( (modes(m,1)*L_x1*x1 + modes(m,2)*L_x2*x2) );
        mm = mm / sqrt(sum((mm(:)) .* conj((mm(:))) )*(L1*L2)/N^2);
        for n = 1:size(modes,1)
            nn = cos( (modes(n,1)*L_x1*x1 + modes(n,2)*L_x2*x2) );
            nn = nn / sqrt(sum((nn(:)) .* conj((nn(:))) )*(L1*L2)/N^2);
            orthocheck(m,n) = sum((mm(:)) .* conj((nn(:))) )*(L1*L2)/N^2;
        end
    end

    switch type
        case 'active'
            for k = 1:size(modes,1)
                if abs(ampstars(k,1)) < 1e-3
                    ampstars(k,1) = NaN;
                    modes(k,:) = [ NaN NaN ];
                end
            end
    end
    modes = rmmissing(modes);
    ampstars = rmmissing(ampstars);

    u_eff = 0;
    for k = 1:size(modes,1)
        u_ef = cos( (modes(k,1)*L_x1*x1 + modes(k,2)*L_x2*x2) );
        u_ef = u_ef / sqrt(sum((u_ef(:)) .* conj((u_ef(:))) )*(L1*L2)/N^2);
        u_ef = ampstars(k,1)*u_ef;
        %[~, max_idx] = max(u_ef(:));
        [~, max_idx] = min(u_ef(:));
        [rowmax, colmax] = ind2sub(size(u_ef), max_idx);
        u1_ef = [ u_ef(rowmax+1:end,:) ;  u_ef(1:rowmax,:)  ];
        u2_ef = [ u1_ef(:,colmax:end) ,  u1_ef(:,1:colmax-1)  ];

        u_eff = u_eff + u2_ef;
        %u_eff = u_eff + u_ef;
    end

    u_centered = u2 - mean(u2(:));
    u_norm = u_centered / sqrt( sum((u_centered(:)) .* conj((u_centered(:))))*(L1*L2)/N^2 );
    phi_norm = u_eff / sqrt( sum((u_eff(:)) .* conj((u_eff(:))))*(L1*L2)/N^2 );
    match_score = 1 - sum((u_norm(:)) .* conj((phi_norm(:))))*(L1*L2)/N^2 ;
    %fprintf('Correlation match: %.7f\n', match_score(1)); 
    
    %{
    h = figure;
    set(gcf,'Position',[100 100 1800 750])
    set(gcf,'color','white')
    set(gca,'color','white') 
    subplot(1,2,1)
    pcolor(x1p, x2p, u2); 
    shading interp;
    xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex','FontSize',18); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex','FontSize',18);
    title(['Phase-shifted $\widetilde{\varphi}_{K,2\pi(' num2str(L_s1) '),' num2str(T) '}$'],'Interpreter','latex','FontSize',18);
    colormap(redblue)
    
    subplot(1,2,2)
    pcolor(x1p, x2p, u_eff);
    shading interp;
    xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex','FontSize',18); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex','FontSize',18);
    switch type
        case {'active','full'}
            title(['$\phi_*(' num2str(L_s1) ',' num2str(L_s2) ')$ with ' num2str(length(ampstars)) ' modes'],'Interpreter','latex','FontSize',18);
        case 'dominant'
            title(['$\phi^d_*(' num2str(L_s1) ',' num2str(L_s2) ')$ with ' num2str(length(ampstars)) ' modes'],'Interpreter','latex','FontSize',18);
    end    
    colormap(redblue)

    filename = [ pwd '/match' num2str(L_s1) '_' IC '_T' num2str(T) '_N' num2str(N)];
    switch type
        case 'active'
            filename = [ pwd '/matcha' num2str(L_s1) '_' IC '_T' num2str(T) '_N' num2str(N)];
        case 'dominant'
            filename = [ pwd '/matchd' num2str(L_s1) '_' IC '_T' num2str(T) '_N' num2str(N)];
    end
    writematrix(match_score, [filename '_matchscore.dat']);
    writematrix(modes, [filename '_modes.dat']);
    writematrix(ampstars, [filename '_amps.dat']);
    saveas(h,[filename '.fig'])
    exportgraphics(h,[filename '.pdf'])
    %}
end
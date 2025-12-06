function process_gif(IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
    energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
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
    fig = figure('Visible', 'on');
    set(fig, 'Position', [100 100 1600 900], 'Color', 'white', 'Resize', 'off');
    
    set(groot, 'DefaultAxesLooseInset', [0 0 0 0]);
    
    ax(1) = subplot(2,4,1);
    ax(2) = subplot(2,4,2);
    ax(3) = subplot(2,4,3);
    ax(4) = subplot(2,4,4);
    ax(5) = subplot(2,4,5);
    ax(6) = subplot(2,4,6);
    ax(7) = subplot(2,4,7);
    
    for k = 1:7
        set(ax(k), ...
            'Color','white', ...
            'FontSize',16, ...
            'LabelFontSizeMultiplier', 1, ...
            'TitleFontSizeMultiplier', 1);
        axis(ax(k), 'square');
        pos = get(ax(k),'Position');
        pos(2) = pos(2) - 0.05;
        set(ax(k),'Position',pos);
    end
    
    ymin_energy = 0.5 * min(energyL2);
    ymax_energy = 1.5 * max(energyH2);
    
    ymax_spec = 1.5 * max(v_mean(:));
    ymax_proj = 1.5 * max(projcoeffradialevolution(:));
    
    title1 = 'Forward-time 2DKS solution';
    title2 = '';
    sg = sgtitle({title1, title2}, 'Interpreter','latex','FontSize',22);
    
    Ntime_remaining = Ntime;
    if utility2(3) > 0
        frames = utility2(3);
    else
        frames = Ntime;
    end
    frameovermax = 0;
    
    switch IC
        case {'optimized'}
            frames = 100;
    end
    
    filename = [pwd '/media/movies/phys_' parameterlist '_frames_' num2str(frames) '.gif'];
    switch IC
        case {'optimized'}
            filename = [pwd '/media/movies/phys_' optparameters '_frames_' num2str(frames) '.gif'];
    end
    
    drow = 0;
    dcol = 0;
    
    % ---------- NEW: handles for graphics objects ----------
    h_surf1 = []; h_surf2 = [];
    h_H2 = []; h_H1 = []; h_L2 = []; h_xline = [];
    a_strip = []; a_xline = []; a_resline = [];
    h_spec1 = []; h_proj = []; h_spec2 = [];
    proj_strip = []; proj_coeffpts = [];
    didInitPlots = false;
    % ------------------------------------------------------
    
    
    for i = 1 : ceil(Ntime/frames) : Ntime+1
    
        if i == Ntime + 1
            i = Ntime;
        end
    
        % ---------- loading logic unchanged ----------
        if Ntime < Ntime_save_max && i == 1
            [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime T], utility1);
        else
            if (Ntime_remaining >= Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                frameovermax = frameovermax + 1;
                currentT = frameovermax*Ntime_save_max/Ntime*T;
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, [Ntime_save_max T], utility1);
                Ntime_remaining = Ntime_remaining - Ntime_save_max;
            elseif (Ntime_remaining < Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                frameovermax = frameovermax + 1;
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, [Ntime_remaining T], utility1);
            end
        end
        % --------------------------------------------
    
        currentT = (i-1)/(Ntime-1)*T;
    
        if mod(i,Ntime_save_max) ~= 0
            imod = mod(i,Ntime_save_max);
        else
            imod = Ntime_save_max;
        end
    
        u_i = reshape(u_n(:,imod), [N, N]);
    
        if i == 1
            [~, min_idx] = min(u_i(:));
            [rowmin, colmin] = ind2sub(size(u_i), min_idx);
            [Ny, Nx] = size(u_i);
            rowc = floor((Ny+1)/2) + 1;
            colc = floor((Nx+1)/2) + 1;
            drow = rowc - rowmin;
            dcol = colc - colmin;
        end
        u_i_ps = circshift(u_i, [drow, dcol]);
    
        u_i2x = [u_i_ps , u_i_ps ; u_i_ps , u_i_ps];
    
        % ==========================================================
        %   FIRST FRAME: CREATE PLOTS & LABELS, STORE HANDLES
        %   LATER FRAMES: ONLY UPDATE DATA, NO cla, NO new plots
        % ==========================================================
        if ~didInitPlots
            % ------------ AXIS 1: physical field ---------------
            h_surf1 = surfc(ax(1), x1pi, x2pi, u_i_ps);
            xlabel(ax(1),'$x_1$','Interpreter','latex','FontSize',16);
            ylabel(ax(1),'$x_2$','Interpreter','latex','FontSize',16);
            shading(ax(1),'interp');
            colormap(ax(1), redblue);
            view(ax(1),3);
            axis(ax(1),'square');
    
            % ------------ AXIS 2: tiled domain -----------------
            h_surf2 = surfc(ax(2), x12x, x22x, u_i2x);
            xlabel(ax(2),'$\frac{x_1}{2\pi}$','Interpreter','latex','FontSize',16);
            ylabel(ax(2),'$\frac{x_2}{2\pi}$','Interpreter','latex','FontSize',16);
            shading(ax(2),'interp');
            colormap(ax(2), redblue);
            xline(ax(2), L_s1, '--');
            yline(ax(2), L_s2, '--');
            view(ax(2),2);
            axis(ax(2),'square');
    
            % ------------ AXIS 3: energy evolution -------------
            h_H2 = semilogy(ax(3), timewindow, energyH2, 'g'); hold(ax(3),'on');
            h_H1 = semilogy(ax(3), timewindow, energyH1, 'r');
            h_L2 = semilogy(ax(3), timewindow, energyL2, 'b');
            h_xline = xline(ax(3), currentT, '-');
            hold(ax(3),'off');
    
            xlabel(ax(3),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(3),'$\| \phi(t;\varphi) \|^2_{S}$','Interpreter','latex','FontSize',16);
            xlim(ax(3), [0 T]);
            ylim(ax(3), [ymin_energy ymax_energy]);
            title(ax(3), "Energy evolution", 'Interpreter','latex','FontSize',16);
            legend(ax(3), '$S=H^2$','$S=H^1$','$S=L^2$', ...
                   'Interpreter','latex','Location','southeast','FontSize',12);
            axis(ax(3),'square');
    
            % ------------ AXIS 4: radial spectrum --------------
            h_spec1 = semilogy(ax(4), v_mean(:,i), "o--");
            xlabel(ax(4),'$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',16); 
            ylabel(ax(4),'$E(k)$','Interpreter','latex','FontSize',16);
            title(ax(4),"Energy spectrum",'Interpreter','latex','FontSize',16);
            xlim(ax(4), [1 size(v_mean,1)]);
            ylim(ax(4), [1e-20 ymax_spec]);
            axis(ax(4),'square');
            
            % ------------ AXIS 5: analyticity strip -----------
            a_strip = plot(ax(5), timewindow, astripwidth, 'b'); 
            hold(ax(5),'on');
            a_xline = xline(ax(5), currentT, '-');
            a_resline = yline(ax(5), max(L1/N,L2/N), '--');
            ylim(ax(5), [0 1.5*max(astripwidth)]);
            hold(ax(5),'off');
    
            xlabel(ax(5),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(5),'$\delta(t)$','Interpreter','latex','FontSize',16);
            title(ax(5), "Analyticity strip width", 'Interpreter','latex','FontSize',16);
            axis(ax(5),'square');

            % ------------ AXIS 6: projection coeffs v radius ------------
            idx0 = find(projcoeffradialevolution(:,end) == 0, 3, 'first');
            idx0 = idx0(end);
            h_proj = plot(ax(6), projcoeffradialevolution(1:idx0,i), "o--");
            xlabel(ax(6),'$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',16); 
            ylabel(ax(6),'$\sum_{k} |a_k|$','Interpreter','latex','FontSize',16);
            title(ax(6),"Projection coefficient spectrum",'Interpreter','latex','FontSize',16);
            xlim(ax(6), [1 idx0]);
            ylim(ax(6), [1e-20 ymax_proj]);
            axis(ax(6),'square');

            % ------------ AXIS 7: projection coeffs v time ------------
            proj_strip = plot(ax(7), timewindow, projcoeffmodeevolution(:,3:end)'); 
            hold(ax(7),'on');
            proj_coeffpts = plot(ax(7), currentT, projcoeffmodeevolution(:,i+2)', 'o');
            ylim(ax(7), [0.99*min(min(projcoeffmodeevolution(:,3:end))) 1.01*max(max(projcoeffmodeevolution(:,3:end)))]);
            hold(ax(7),'off');
    
            xlabel(ax(7),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(7),'$P\left(a_k\right)$','Interpreter','latex','FontSize',16);
            %ylabel(ax(7),'${|a_k|}/{\| \sum |a_k| \|}$','Interpreter','latex','FontSize',16);
            title(ax(7), "Projection coefficient evolution", 'Interpreter','latex','FontSize',16);
            axis(ax(7),'square');
    
            didInitPlots = true;
    
        else
            % ================= UPDATE ONLY =====================
            % Axis 1 & 2 surfaces:
            set(h_surf1, 'ZData', u_i_ps);
            set(h_surf2, 'ZData', u_i2x);
    
            % Axis 3: just move the xline (curves are static arrays)
            h_xline.Value = currentT;
    
            % Axis 4: update spectrum
            set(h_spec1, 'YData', v_mean(:,i));
                    
            % Axis 5: update analyticity strip
            a_xline.Value = currentT;

            % Axis 6: update projection coefficient spectrum
            set(h_proj, 'YData', projcoeffradialevolution(:,i));

            % Axis 7: update proj coeff evolution
            set(proj_coeffpts, 'XData', currentT*ones(size(projcoeffmodeevolution,1),1), 'YData', projcoeffmodeevolution(:,i+2));
    
            % ==================================================
        end
    
        % ------------ Update LaTeX sgtitle --------------------
        if strcmp(IC,'optimized')
            phi_str = '\tilde{\varphi}';
        else
            phi_str = ['\varphi_{' IC '}'];
        end
    
        title2 = sprintf( ...
            '$\\varphi = %s, N = %d, {\\Delta}t = %.5g, K = %.0f, L_1 = 2\\pi(%.2f), L_2 = 2\\pi(%.2f), T = %.5f$', ...
            phi_str, N, dt, K, L_s1, L_s2, currentT);
    
        sg.String = {title1, title2};
    
        drawnow limitrate;   % lighter than bare drawnow
    
        if i == 1
            gif(filename,'overwrite',true);
        else
            gif;
        end
    
    end

end
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
    projcoeffradialcut = projcoeffradialevolution(~any(projcoeffradialevolution == 0, 2), :);
    num_radii = size(projcoeffradialcut,1);
    number_of_plots = 7 + num_radii;
    maxplotcols = 5;
    number_of_plot_rows = ceil(number_of_plots/maxplotcols);
    figwidth = 400*maxplotcols;
    figheight = 450*number_of_plot_rows;

    fig = figure('Visible', 'on');
    set(fig, 'Position', [100 100 figwidth figheight], 'Color', 'white', 'Resize', 'off');
    
    set(groot, 'DefaultAxesLooseInset', [0 0 0 0]);
    
    for k = 1:number_of_plots
        ax(k) = subplot(number_of_plot_rows,maxplotcols,k);
    end
    
    for k = 1:number_of_plots
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
            k = 1;
            h_surf1 = surfc(ax(k), x1pi, x2pi, u_i_ps);
            xlabel(ax(k),'$x_1$','Interpreter','latex','FontSize',16);
            ylabel(ax(k),'$x_2$','Interpreter','latex','FontSize',16);
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            view(ax(k),3);
            axis(ax(k),'square');
    
            % ------------ AXIS 2: tiled domain -----------------
            k = 2;
            h_surf2 = surfc(ax(k), x12x, x22x, u_i2x);
            xlabel(ax(k),'$\frac{x_1}{2\pi}$','Interpreter','latex','FontSize',16);
            ylabel(ax(k),'$\frac{x_2}{2\pi}$','Interpreter','latex','FontSize',16);
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            xline(ax(k), L_s1, '--');
            yline(ax(k), L_s2, '--');
            view(ax(k),2);
            axis(ax(k),'square');
    
            % ------------ AXIS 3: energy evolution -------------
            k = 3;
            h_H2 = semilogy(ax(k), timewindow, energyH2, 'g'); hold(ax(3),'on');
            h_H1 = semilogy(ax(k), timewindow, energyH1, 'r');
            h_L2 = semilogy(ax(k), timewindow, energyL2, 'b');
            %h_xline = xline(ax(k), currentT, '-');
            H2_xline = plot(ax(k), timewindow(1,i), energyH2(1,i), 'ko');
            H1_xline = plot(ax(k), timewindow(1,i), energyH1(1,i), 'ko');
            L2_xline = plot(ax(k), timewindow(1,i), energyL2(1,i), 'ko');
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(k),'$\| \phi(t;\varphi) \|^2_{S}$','Interpreter','latex','FontSize',16);
            xlim(ax(k), [0 T]);
            ylim(ax(k), [ymin_energy ymax_energy]);
            title(ax(k), "Energy evolution", 'Interpreter','latex','FontSize',16);
            legend(ax(k), '$S=H^2$','$S=H^1$','$S=L^2$', ...
                   'Interpreter','latex','Location','southeast','FontSize',12);
            axis(ax(k),'square');
    
            % ------------ AXIS 4: radial spectrum --------------
            k = 4;
            h_spec1 = semilogy(ax(k), v_mean(:,1), v_mean(:,i+1), "-");
            xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',16); 
            ylabel(ax(k),'$E(k)$','Interpreter','latex','FontSize',16);
            title(ax(k),"Energy spectrum",'Interpreter','latex','FontSize',16);
            xlim(ax(k), [1 size(v_mean,1)]);
            ylim(ax(k), [1e-20 ymax_spec]);
            axis(ax(k),'square');
            
            % ------------ AXIS 5: analyticity strip -----------
            k = 5;
            a_strip = plot(ax(k), timewindow, astripwidth, 'b'); 
            hold(ax(k),'on');
            %a_xline = xline(ax(k), currentT, '-');
            a_xline = plot(ax(k), timewindow(1,i), astripwidth(1,i), 'ko');
            a_resline = yline(ax(k), max(L1/N,L2/N), '--');
            ylim(ax(k), [0 1.5*max(astripwidth)]);
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(k),'$\delta(t)$','Interpreter','latex','FontSize',16);
            title(ax(k), "Analyticity strip width", 'Interpreter','latex','FontSize',16);
            axis(ax(k),'square');

            % ------------ AXIS 6: projection coeffs v radius ------------
            k = 6;
            %{
            idx0 = find(projcoeffradialevolution(:,end) == 0, 1, 'first');
            idx0 = idx0(end);
            projcoeffradialevolution(idx0,1) = round(projcoeffradialevolution(idx0-1,1))+1;
            h_proj = plot(ax(k), projcoeffradialevolution(1:idx0,1), projcoeffradialevolution(1:idx0,i+1), "o--");
            %}
            h_proj = plot(ax(k), projcoeffradialcut(:,1), projcoeffradialcut(:,i+1), "o--");
            xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',16); 
            ylabel(ax(k),'$P(a_k)$','Interpreter','latex','FontSize',16);
            title(ax(k),"Projection coefficient spectrum",'Interpreter','latex','FontSize',16);
            %xlim(ax(k), [1 idx0]);
            ylim(ax(k), [min(min(projcoeffradialcut(:,2:end))) max(max(projcoeffradialcut(:,2:end)))]);
            axis(ax(k),'square');

            % ------------ AXIS 7-X: projection coeffs v time ------------
            k = 7;
            proj_strip(k,1:size(projcoeffmodeevolution,1)) = plot(ax(k), timewindow, projcoeffmodeevolution(:,3:end)'); 
            hold(ax(k),'on');
            
            Nmodes = size(projcoeffmodeevolution,1);
            % Initial coordinates (for the first frame)
            xpts = currentT * ones(Nmodes,1);
            ypts = projcoeffmodeevolution(:, i+2);
            %kx = projcoeffmodeevolution(:,1);
            %ky = projcoeffmodeevolution(:,2);
            %labels = compose('%g', 1:Nmodes);%compose('%g , %g', kx, ky);

            % Create marker handles (one per point)
            proj_coeffpts(k,1:Nmodes) = plot(ax(k), xpts, ypts, 'ko');
            hold(ax(k), 'on')
            
            % Create text labels (one per point)
            %{
            for r = 1:Nmodes
                proj_labels(k,r) = text(ax(k), xpts(r), ypts(r), labels(r), ...
                    'VerticalAlignment',  'middle', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize',            8);
            end
            %}

            ylim(ax(k), [0.99*min(min(projcoeffmodeevolution(:,3:end))) 1.01*max(max(projcoeffmodeevolution(:,3:end)))]);
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$','Interpreter','latex','FontSize',16);
            ylabel(ax(k),'$P\left(a_k\right)$','Interpreter','latex','FontSize',16);
            %ylabel(ax(k),'${|a_k|}/{\| \sum |a_k| \|}$','Interpreter','latex','FontSize',16);
            title(ax(k), 'Full projection', 'Interpreter','latex','FontSize',16);
            legend(ax(k), proj_coeffpts(k,1), [num2str(Nmodes) ' modes'],'Location','southwest')
            axis(ax(k),'square');
                        
            radnum = 0;
            averagemodes = size(projcoeffmodeevolution,1)/size(rmmissing(projcoeffradialevolution),1);
            radrows = NaN(num_radii,averagemodes);
            for k = 8:number_of_plots
                radnum = radnum + 1;
                radius = projcoeffradialevolution(radnum,1);
                currows = find( abs(vecnorm(projcoeffmodeevolution(:,1:2), 2, 2) - radius) < 1e-10 );
                radrows(radnum,1:size(currows,1)) = currows;
                proj_strip(k,1:size(currows,1)) = plot(ax(k), timewindow, projcoeffmodeevolution(radrows(radnum,1:size(currows,1)),3:end)'); 
                hold(ax(k),'on');
                proj_coeffpts(k,1:size(currows,1)) = plot(ax(k), currentT*ones(size(radrows(radnum,1:size(currows,1)),1),1), projcoeffmodeevolution(radrows(radnum,1:size(currows,1)),i+2)', 'ko');
                %ylim(ax(k), [0.99*min(min(projcoeffmodeevolution(radrows(radnum,1:size(currows,1)),3:end))) 1.01*max(max(projcoeffmodeevolution(radrows(radnum,1:size(currows,1)),3:end)))]);
                ylim(ax(k), [0.99*min(min(projcoeffmodeevolution(:,3:end))) 1.01*max(max(projcoeffmodeevolution(:,3:end)))]);
                
                hold(ax(k),'off');
        
                xlabel(ax(k),'Time $t$','Interpreter','latex','FontSize',16);
                ylabel(ax(k),'$P\left(a_k\right)$','Interpreter','latex','FontSize',16);
                %ylabel(ax(k),'${|a_k|}/{\| \sum |a_k| \|}$','Interpreter','latex','FontSize',16);
                title(ax(k), ['$k =$' num2str(radius) ' projection'], 'Interpreter','latex','FontSize',16);
                legend(ax(k), proj_coeffpts(k,1), [num2str(size(currows,1)) ' modes'],'Location','southwest')
                axis(ax(k),'square');
            end
    
            didInitPlots = true;
    
        else
            % ================= UPDATE ONLY =====================
            % Axis 1 & 2 surfaces:
            set(h_surf1, 'ZData', u_i_ps);
            set(h_surf2, 'ZData', u_i2x);
    
            % Axis 3: just move the xline (curves are static arrays)
            %h_xline.Value = currentT;
            set(H2_xline, 'XData', timewindow(1,i), 'YData', energyH2(i,1));
            set(H1_xline, 'XData', timewindow(1,i), 'YData', energyH1(i,1));
            set(L2_xline, 'XData', timewindow(1,i), 'YData', energyL2(i,1));
    
            % Axis 4: update spectrum
            set(h_spec1, 'YData', v_mean(:,i+1));
                    
            % Axis 5: update analyticity strip
            %a_xline.Value = currentT;
            set(a_xline, 'XData', timewindow(1,i), 'YData', astripwidth(i,1));

            % Axis 6: update projection coefficient spectrum
            %set(h_proj, 'YData', projcoeffradialevolution(1:idx0,i+1));
            set(h_proj, 'YData', projcoeffradialcut(:,i+1));

            % Axis X: update proj coeff evolution
            k = 7;
            xpts = currentT * ones(Nmodes,1);
            ypts = projcoeffmodeevolution(:, i+2);
            set(proj_coeffpts(k,1:Nmodes), 'XData', xpts, 'YData', ypts);

            radnum = 0;
            for k = 8:number_of_plots
                radnum = radnum + 1;
                numrows = size(rmmissing(radrows(radnum,:)),2);
                set(proj_coeffpts(k,1:numrows), 'XData', currentT*ones(numrows,1), 'YData', projcoeffmodeevolution(radrows(radnum,1:numrows),i+2));
            end

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
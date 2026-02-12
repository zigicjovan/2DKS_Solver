function process_gif(u_IC, IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
    energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
    parameterlist,optparameters)

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

    % reshape
    k1vecN = k1_n(:);
    k2vecN = k2_n(:);
    Klap = k1_l.^2 + k2_l.^2;                                   % for 2nd order linear term (Laplacian) 
    kvecl2 = Klap(:);
    KK = (k1_l.^2 + k2_l.^2).^2;                                % for 4th order linear term (bi-Laplacian)
    kvecl4 = KK(:);         
    
    % Fourier space operators
    Lin = 1i^2 * (kvecl2) + 1i^4 * (kvecl4);                    % Linear operator
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

    %% plot grid
    % 4-period physical space domain
    x12_pts = 2*L_s1*linspace( 0 , 1 - 1/N , 2*N ); 
    x22_pts = 2*L_s2*linspace( 0 , 1 - 1/N , 2*N ); 
    [ x12x , x22x ] = meshgrid(x12_pts,x22_pts); % 2-dimensional grid
    %%

    % AS strip width fit
    asstrip_fit = v_mean(:,2:end);
    for i = 2:size(v_mean,2)
        asstrip_fit(:,i-1) = astripwidth(i-1,1)*exp(-2*astripwidth(i-1,2)*v_mean(:,1));
    end

    %% process
    projcoeffradialcut = projcoeffradialevolution(~any(projcoeffradialevolution == 0, 2), :);
    radialcutmodes = NaN(size(projcoeffradialcut,1),2);
    modecounter = 1;
    k1 = 1;
    k2 = 0;
    while modecounter < (size(projcoeffradialcut,1) + 1) 
        currentradius = sqrt(k1^2 + k2^2);
        idx = find( abs( currentradius - projcoeffradialcut(:,1) ) < 1e-10 );
        if ~isempty(idx) && isnan(radialcutmodes(idx,1))
            radialcutmodes(idx,:) = [ k1 , k2 ];
            modecounter = modecounter + 1;
            if currentradius <= projcoeffradialcut(end,1)
                k1 = k1 + 1;
            else
                k1 = k2 + 1;
                k2 = k2 + 1;
            end
        else
            if currentradius <= projcoeffradialcut(end,1)
                k1 = k1 + 1;
            else
                k1 = k2 + 1;
                k2 = k2 + 1;
            end
        end
    end
    modelabels = "(" + string(radialcutmodes(:,1)) + "," + string(radialcutmodes(:,2)) + ")";
    modecats = categorical(modelabels, modelabels, 'Ordinal', true);
    %num_radii = size(projcoeffradialcut,1);
    %number_of_plots = 7 + num_radii;
    number_of_plots = 7;
    maxplotcols = 4;
    number_of_plot_rows = ceil(number_of_plots/maxplotcols);
    figwidth = 400*maxplotcols;
    figheight = 450*number_of_plot_rows;

    fig = figure('Visible', 'off');
    set(fig, 'Position', [100 100 figwidth figheight], 'Color', 'white', 'Resize', 'off');
    
    wordsize = 16;

    set(groot, ...
    'DefaultAxesLooseInset',            [0 0 0 0], ...
    'DefaultLegendFontSize',            ceil(3*wordsize/4), ...
    'DefaultTextInterpreter',           'latex', ...
    'DefaultAxesTickLabelInterpreter',  'latex', ...
    'DefaultLegendInterpreter',         'latex');
    
    for k = 1:number_of_plots
        ax(k) = subplot(number_of_plot_rows,maxplotcols,k);
    end
    
    for k = 1:number_of_plots
        set(ax(k),'Color','white');
        axis(ax(k), 'square');
        pos = get(ax(k),'Position');
        pos(2) = pos(2) - 0.05;
        set(ax(k),'Position',pos);
    end
    
    ymin_energy = 0.5 * min(energyL2);
    ymax_energy = 1.5 * max(energyH2);
    
    ymax_spec = 1.5 * max(v_mean(:));
    
    title1 = 'Forward-time 2DKS solution';
    title2 = '';
    sg = sgtitle({title1, title2}, 'Interpreter','latex','FontSize',ceil(1.25*wordsize));
    
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
    h_surf1 = []; h_surf2 = []; a_xline = []; 
    h_spec1 = []; h_proj = []; h_spec2 = [];
    proj_strip = []; proj_coeffpts = [];
    didInitPlots = false;
    % ------------------------------------------------------
    for i = 1:Ntime

        % nonlinear terms and solution substeps
        if i > 1
            v_1 = ks_nonlinear_fwd( v_step, D1vec, D2vec, N, Lin, dt, alpha_I, beta_I, alpha_E, beta_E );
        else
            v_1 = v_step;
        end

        % correction of non-zero mean solution
        if abs(mean(v_1)) > 1e-5 
            v_1(1) = 0;
        end

        % save solution step to workspace
        v_12 = reshape( v_1, [ N , N ] );
        u_n = real(ifft2(v_12));  
        v_step = v_1;

        currentT = (i-1)/(Ntime-1)*T;
        u_i = u_n;
    %{
    for i = 1 : ceil(Ntime/frames) : Ntime+1
    
        if i == Ntime + 1
            i = Ntime;
        end
        % ---------- reading from file ----------
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

        if mod(i,Ntime_save_max) ~= 0
            imod = mod(i,Ntime_save_max);
        else
            imod = Ntime_save_max;
        end
    
        u_i = reshape(u_n(:,imod), [N, N]);
        %}
    
        % appropriate phase-shift to center any local minima
        if i == 1
            minctr = 1;
            [~,minidx] = mink(u_i(:),20);
            while minctr < size(minidx,1)
                minidx = [minidx(1:minctr) ; minidx(abs(minidx(minctr) - minidx(minctr:end)) > .05*size(u_i(:),1))];
                minctr = minctr + 1;
            end
            [row, col] = ind2sub(size(u_i), minidx);
            sz = size(u_i);
            
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
            h_surf1 = surfc(ax(k), x1, x2, u_i_ps);
            xlabel(ax(k),'$x_1$');
            ylabel(ax(k),'$x_2$');
            %title(ax(k),"Periodic solution field");
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            view(ax(k),3);
            axis(ax(k),'square');
    
            % ------------ AXIS 2: tiled domain -----------------
            k = 2;
            h_surf2 = surfc(ax(k), x12x, x22x, u_i2x);
            xlabel(ax(k),'$\frac{x_1}{2\pi}$' );
            ylabel(ax(k),'$\frac{x_2}{2\pi}$' );
            title(ax(k),"Tiled solution field" );
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            xline(ax(k), L_s1, '--');
            yline(ax(k), L_s2, '--');
            view(ax(k),2);
            axis(ax(k), [0 2*L_s1 0 2*L_s2])
            axis(ax(k), 'equal')
    
            % ------------ AXIS 3: energy evolution -------------
            k = 3;
            semilogy(ax(k), timewindow, energyH2, 'g'); 
            hold(ax(k),'on');
            semilogy(ax(k), timewindow, energyH1, 'r');
            semilogy(ax(k), timewindow, energyL2, 'b');
            %h_xline = xline(ax(k), currentT, '-');
            H2_xline = plot(ax(k), timewindow(1,i), energyH2(i,1), 'ko');
            H1_xline = plot(ax(k), timewindow(1,i), energyH1(i,1), 'ko');
            L2_xline = plot(ax(k), timewindow(1,i), energyL2(i,1), 'ko');
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$' );
            ylabel(ax(k),'$\| \phi(t;\varphi) \|^2_{S}$' );
            xlim(ax(k), [0 T]);
            ylim(ax(k), [ymin_energy ymax_energy]);
            title(ax(k), "Energy evolution");
            legend(ax(k), '$S=H^2$','$S=H^1$','$S=L^2$','Location','southeast','Box','off');
            axis(ax(k),'square');
    
            % ------------ AXIS 4: radial spectrum --------------
            k = 4;
            h_spec1 = semilogy(ax(k), v_mean(:,1), v_mean(:,i+1), "-");
            hold(ax(k),'on');
            h_spec2 = semilogy(ax(k), v_mean(:,1), asstrip_fit(:,i), "r--");
            hold(ax(k),'off');
            xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$' ); 
            legend(ax(k),'$E(k)$','$Ce^{-2\delta k}$','Box','off','FontSize',wordsize);
            title(ax(k),"Energy spectrum" );
            xlim(ax(k), [1 size(v_mean,1)]);
            ylim(ax(k), [1e-20 ymax_spec]);
            axis(ax(k),'square');
            
            % ------------ AXIS 5: analyticity strip -----------
            k = 5;
            plot(ax(k), timewindow, astripwidth(:,2), 'b'); 
            hold(ax(k),'on');
            %a_xline = xline(ax(k), currentT, '-');
            a_xline = plot(ax(k), timewindow(1,i), astripwidth(i,2), 'ko');
            yline(ax(k), max(L1/N,L2/N), '--');
            ylim(ax(k), [0 1.5*max(astripwidth(:,2))]);
            xlim(ax(k), [0 T]);
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$' );
            ylabel(ax(k),'$\delta(t)$' );
            title(ax(k), "Analyticity strip width");
            axis(ax(k),'square');

            % ------------ AXIS 6: projection coeffs v radius ------------
            k = 6;
            %{
            idx0 = find(projcoeffradialevolution(:,end) == 0, 1, 'first');
            idx0 = idx0(end);
            projcoeffradialevolution(idx0,1) = round(projcoeffradialevolution(idx0-1,1))+1;
            h_proj = plot(ax(k), projcoeffradialevolution(1:idx0,1), projcoeffradialevolution(1:idx0,i+1), "o--");
            %}
            cla(ax(k)) 
            xmodelabels = 1:numel(modecats);
            h_proj = bar(ax(k), xmodelabels, projcoeffradialcut(:, i+1), 'BarWidth', 1.0);
            if numel(modecats) < 15
                ax(k).XTick = xmodelabels;
                ax(k).XTickLabel = cellstr(modecats);
                xlabel(ax(k),'$(k_1,k_2)$' );
            else
                xlabel(ax(k),'$(k_1,0)$' );
            end
            ax(k).XLim = [0.5, numel(modecats) + 0.5];
            ylabel(ax(k),'$P(a_k)$' );
            title(ax(k),"Projection coefficient weights" );
            %xlim(ax(k), [1 idx0]);
            ylim(ax(k), [0.9*min(min(projcoeffradialcut(:,2:end))) , 1.1*max(max(projcoeffradialcut(:,2:end)))]);
            axis(ax(k),'square');

            % ------------ AXIS 7-X: projection coeffs v time ------------
            k = 7;
            proj_strip(k,1:size(projcoeffmodeevolution,1)) = plot(ax(k), timewindow, projcoeffmodeevolution(:,3:end)'); 
            modalfigylim = [0.99*min(min(projcoeffmodeevolution(:,3:end))) 1.01*max(max(projcoeffmodeevolution(:,3:end)))];
            hold(ax(k),'on');
            
            Nmodes = size(projcoeffmodeevolution,1);
            % Initial coordinates (for the first frame)
            projcoeff_xpts = currentT * ones(Nmodes,1);
            projcoeff_ypts = [ sqrt(projcoeffmodeevolution(:,1).^2 + projcoeffmodeevolution(:,2).^2) , projcoeffmodeevolution(:, i+2) ];
            %kx = projcoeffmodeevolution(:,1);
            %ky = projcoeffmodeevolution(:,2);
            %labels = compose('%g', 1:Nmodes);%compose('%g , %g', kx, ky);

            % Create marker handles (one per point)
            proj_coeffpts(k,1:Nmodes) = plot(ax(k), projcoeff_xpts, projcoeff_ypts(:,2), 'ko');
            hold(ax(k), 'on')
            
            % Create text labels (one per radius) and order to prevent overlap
            modelabelposition = cell(size(modelabels,1),3);
            figfrac = (modalfigylim(2)-modalfigylim(1))*0.07;
            for imode = 1:size(modelabels,1)
                idx = find(( abs( projcoeff_ypts(:,1) - projcoeffradialcut(imode,1) ) < 1e-10 ) , 1 , 'first');
                modelabelposition{imode,1} = modelabels(imode);
                modelabelposition{imode,2} = idx;
                modelabelposition{imode,3} = projcoeffmodeevolution(idx,end);
            end
            modelabelposition = sortrows(modelabelposition,3);
            if size(modelabels,1) > 15
                modelabelposition{1,3} = modalfigylim(1);
            end
            for imode = 2:size(modelabels,1)
                posdiff = modelabelposition{imode,3} - modelabelposition{imode-1,3};
                if posdiff < figfrac && size(modelabels,1) < 15
                    modelabelposition{imode,3} = modelabelposition{imode,3} + figfrac - posdiff;
                elseif size(modelabels,1) > 15
                    modelabelposition{imode,3} = modelabelposition{imode-1,3} + figfrac;
                end
            end
            if numel(modecats) < 15
                for imode = 1:size(modelabels,1)
                    text(ax(k), T, modelabelposition{imode,3}, modelabelposition{imode,1},...
                        'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                        'Interpreter','latex', 'FontSize', ceil(0.75*wordsize));
                end
            end

            xlim(ax(k), [0 T]);
            ylim(ax(k), modalfigylim);
            hold(ax(k),'off');
    
            xlabel(ax(k),'Time $t$' );
            ylabel(ax(k),'$P\left(a_k\right)$' );
            title(ax(k), 'Modal Energy');
            %legend(ax(k), proj_coeffpts(k,1), [num2str(Nmodes) ' modes'],'Location','southwest')
            axis(ax(k),'square');
                     
            %{
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
        
                xlabel(ax(k),'Time $t$' );
                ylabel(ax(k),'$P\left(a_k\right)$' );
                %ylabel(ax(k),'${|a_k|}/{\| \sum |a_k| \|}$' );
                title(ax(k), ['$k =$' num2str(radius) ' projection']);
                legend(ax(k), proj_coeffpts(k,1), [num2str(size(currows,1)) ' modes'],'Location','southwest')
                axis(ax(k),'square');
            end
            %}
    
            % ---- enforce fonts/interpreters AFTER all plotting calls ----
            for kk = 1:number_of_plots
                set(ax(kk), ...
                    'FontSize', ceil(0.75*wordsize), ...
                    'LabelFontSizeMultiplier', 1.0, ...
                    'TitleFontSizeMultiplier', 1.0, ...
                    'TickLabelInterpreter', 'latex');  % axes tick labels
            
                % Ensure existing labels/titles use LaTeX + correct size
                ax(kk).XLabel.Interpreter = 'latex';
                ax(kk).YLabel.Interpreter = 'latex';
                ax(kk).Title.Interpreter  = 'latex';
            
                ax(kk).XLabel.FontSize = wordsize;
                ax(kk).YLabel.FontSize = wordsize;
                ax(kk).Title.FontSize  = wordsize;
            end
            
            % Legends: force font size/interpreter (since legend does not inherit axes)
            lg = findall(fig, 'Type', 'Legend');
            set(lg, 'Interpreter','latex', 'FontSize', ceil(0.75*wordsize));
            % ------------------------------------------------------------

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

            didInitPlots = true;
    
        elseif mod(i,ceil(Ntime/frames)) == 0
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
            set(h_spec2, 'YData', asstrip_fit(:,i));
                    
            % Axis 5: update analyticity strip
            %a_xline.Value = currentT;
            set(a_xline, 'XData', timewindow(1,i), 'YData', astripwidth(i,2));

            % Axis 6: update projection coefficient spectrum
            %set(h_proj, 'YData', projcoeffradialevolution(1:idx0,i+1));
            set(h_proj, 'YData', projcoeffradialcut(:,i+1));

            % Axis X: update proj coeff evolution
            k = 7;
            projcoeff_xpts = currentT * ones(Nmodes,1);
            projcoeff_ypts(:,2) = projcoeffmodeevolution(:, i+2);
            set(proj_coeffpts(k,1:Nmodes), 'XData', projcoeff_xpts, 'YData', projcoeff_ypts(:,2));

            %{
            radnum = 0;
            for k = 8:number_of_plots
                radnum = radnum + 1;
                numrows = size(rmmissing(radrows(radnum,:)),2);
                set(proj_coeffpts(k,1:numrows), 'XData', currentT*ones(numrows,1), 'YData', projcoeffmodeevolution(radrows(radnum,1:numrows),i+2));
            end
            %}

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
end

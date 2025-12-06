function process_figure(figuretype,originalIC, IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,tol,... 
    energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
    parameterlist,optparameters,parfiglist,optparfiglist)

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

    switch figuretype
        case 'energy'
            % Wavenumber evolution plot
            timewindow = linspace(0,T,Ntime);
            h = figure('Visible', 'off');
            semilogy(timewindow,v_mean(1,:),'LineWidth',0.1,'Marker','.')
            hold on;
            for i = 2:size(v_mean,1)
                semilogy(timewindow,v_mean(i,:),'LineWidth',0.1,'Marker','.')
            end
            set(gcf,'Position',[100 100 900 750])
            xlabel('Time $t$','Interpreter','latex'); 
            xlim([0 T])
            ylim([1e-20 max(v_mean(1,:))+1e5 ])
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of Fourier spectrum','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
            %legend("Fourier band", 'Location','southeast','NumColumns',9,'Interpreter','latex')
            %frame = getframe(h);
            %im = frame2im(frame);
            filename = [pwd '/media/energy/wavenumberevol_' parameterlist ];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/energy/wavenumberevol_' optparameters ];
            end
            %imwrite(im,filename,'png')
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])
    
            % Wavenumber IC plot
            h = figure('Visible', 'off');
            semilogy(v_mean(:,1),"o--")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            xlim([1 size(v_mean,1)])
            ylim([ 1e-20 max(v_mean(1,:))+1e5 ])
            set(gcf,'Position',[100 100 900 750])
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Initial Fourier spectrum','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
            filename = [pwd '/media/energy/spectrumIC_' parameterlist ];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/energy/spectrumIC_' optparameters ];
            end
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])
    
            % Wavenumber TC plot
            h = figure('Visible', 'off');
            semilogy(v_mean(:,end),"o--")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            xlim([1 size(v_mean,1)])
            ylim([ 1e-20 max(v_mean(1,:))+1e5 ])
            set(gcf,'Position',[100 100 900 750])
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Terminal Fourier spectrum','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
            filename = [pwd '/media/energy/spectrumTC_' parameterlist ];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/energy/spectrumTC_' optparameters ];
            end
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])
    
            % L2 energy plot
            h = figure('Visible', 'off');
            semilogy(timewindow,energyL2,'LineWidth',0.5,'Marker','.')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Time $t$','Interpreter','latex'); 
            xlim([0 T])
            ylabel('$\| {\phi(t;\varphi)} \|^2_{L^2}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of $L^2$ energy','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
            filename = [pwd '/media/energy/energyL2_' parameterlist ];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/energy/energyL2_' optparameters ];
                otherwise
                    saveas(h,[filename '.fig'])
                    exportgraphics(h,[filename '.pdf'])
            end
        case 'diagnostics'
                  
            energyL2_og = utility1;
            [diagnostics, linesearchJ] = load_2DKSsolution('optimization', 'optimized', dt, T, N, K, L_s1, L_s2, [tol T], originalIC);

            % L2 energy plot
            h = figure('Visible', 'off');
            semilogy(timewindow,energyL2_og,'LineWidth',0.5,'Marker','.')
            hold on
            semilogy(timewindow,energyL2,'LineWidth',0.5,'Marker','.')
            hold off
            set(gcf,'Position',[100 100 900 750])
            xlabel('Time $t$','Interpreter','latex'); 
            xlim([0 T])
            ylabel('$\| {\phi(t;\varphi)} \|^2_{L^2}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of optimized $L^2$ energy','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            legend(['$\phi(t,\varphi_{' originalIC '})$'],'$\phi(t,\tilde{\varphi})$','Interpreter','latex','Location','northwest')
            filename = [pwd '/media/optimization/energyL2comp_' optparameters];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % J history plot
            plotdata = diagnostics(:,1);
            h = figure('Visible', 'off');
            plot(plotdata,'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\mathcal{J}_T(\varphi^{(n)})$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of objective functional $\mathcal{J}_T(\varphi^{(n)})= \| {\phi^{(n)}(T)} \|^2_{L^2}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/Jhistory_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % J change plot
            plotdata = diagnostics(:,2);
            h = figure('Visible', 'off');
            semilogy(plotdata,'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\frac{\mathcal{J}_T(\varphi^{(n+1)}) - \mathcal{J}_T(\varphi^{(n)})}{\mathcal{J}_T(\varphi^{(n)})}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Relative change in objective functional $\mathcal{J}_T(\varphi^{(n)})= \| {\phi^{(n)}(T)} \|^2_{L^2}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/Jchange_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % stepsize plot
            plotdata = diagnostics(:,3);
            h = figure('Visible', 'off');
            semilogy(abs(plotdata),'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\tau^{(n)}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Optimal step-size $\tau^{(n)}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/stepsize_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % manifold plot
            plotdata = diagnostics(:,4);
            h = figure('Visible', 'off');
            semilogy(plotdata,'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\| \varphi^{(n)} \|^2_{L^2}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of initial condition magnitude $K= \| \varphi^{(n)} \|^2_{L^2}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/initialK_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % gradient magnitude plot
            plotdata = diagnostics(:,6);
            h = figure('Visible', 'off');
            semilogy(plotdata,'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\| \nabla \mathcal{J}_T(\varphi^{(n)}) \|^2_{L^2}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of objective gradient magnitude $\| \nabla \mathcal{J}_T(\varphi^{(n)}) \|^2_{L^2}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/Jgradient_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % momentum magnitude plot
            plotdata = abs(diagnostics(:,7));
            h = figure('Visible', 'off');
            semilogy(plotdata,'r-*')
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\|\beta^{(n)}\|^2_{L^2}$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of momentum magnitude $\|\beta^{(n)}\|^2_{L^2}$','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/momentumsize_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            % line search plot
            h = figure('Visible', 'off');
            semilogy(linesearchJ(:,1))
            hold on
            for iter = 2:(length(plotdata)-1)
                semilogy(linesearchJ(:,iter))
            end
            hold off
            set(gcf,'Position',[100 100 900 750])
            xlabel('Iteration number $n$','Interpreter','latex'); 
            xlim([1 length(plotdata)])
            ylabel('$\mathcal{J}_T(\varphi^{(n)})$','Interpreter','latex');
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white')    
            title('Evolution of Brent''s method objective functional evaluation','Interpreter','latex')
            subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            filename = [pwd '/media/optimization/linesearch_' optparameters ];
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])

            close all
        case 'state'
            u_IC = utility1;
            u_TC = utility2;
            %% initial
            h = figure('Visible', 'off');
            axis tight manual % this ensures that getframe() returns a consistent size
            set(gcf,'Position',[100 100 900 750])
    
            % inspect physical solution 
            u_T = reshape( u_IC , [ N , N ] );
            % phase-shift
            [~, min_idx] = min(u_T(:));
            [rowmin, colmin] = ind2sub(size(u_T), min_idx);
            [Ny, Nx] = size(u_T);
            rowc = floor((Ny+1)/2) + 1;
            colc = floor((Nx+1)/2) + 1;
            drow = rowc - rowmin;
            dcol = colc - colmin;
            u_T_ps = circshift(u_T, [drow, dcol]);
    
            % surface plot
            surfc(x1,x2,u_T_ps); 
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
            shading interp
            %pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white') 
            colormap(redblue)
            view(3);
            title('Initial state','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
    
            % contour plot
            view(2);
            axis square;
            % Save image
            filename = [pwd '/media/figures/state/phys_' parameterlist '_initial_contour'];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/figures/state/phys_' optparameters '_initial_contour'];
            end
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])
    
            %% terminal
            h = figure('Visible', 'off');
            axis tight manual % this ensures that getframe() returns a consistent size
            set(gcf,'Position',[100 100 900 750])
    
            % inspect physical solution
            u_T = reshape( u_TC , [ N , N ] );
            % phase-shift
            [~, min_idx] = min(u_T(:));
            [rowmin, colmin] = ind2sub(size(u_T), min_idx);
            [Ny, Nx] = size(u_T);
            rowc = floor((Ny+1)/2) + 1;
            colc = floor((Nx+1)/2) + 1;
            drow = rowc - rowmin;
            dcol = colc - colmin;
            u_T_ps = circshift(u_T, [drow, dcol]);
    
            % surface plot
            surfc(x1,x2,u_T_ps);  
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
            shading interp
            %pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
            view(3);
            fontsize(12,"points")
            set(gca,'fontsize', 16) 
            set(gcf,'color','white')
            set(gca,'color','white') 
            colormap(redblue)
            title('Terminal state','Interpreter','latex')
            subtitle(parfiglist,'Interpreter','latex','FontSize',14)
            switch IC 
                case {'optimized'}
                    subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
            end
    
            % contour plot
            view(2);
            axis square;
            % Save image
            filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal_contour'];
            switch IC
                case {'optimized'}
                    filename = [pwd '/media/figures/state/phys_' optparameters '_terminal_contour'];
            end
            saveas(h,[filename '.fig'])
            exportgraphics(h,[filename '.pdf'])
    end
end
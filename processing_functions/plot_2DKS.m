function maxL2inT = plot_2DKS(save_each, solplot, IC, N, dt, T, K, L_s1, L_s2, Ntime_save_max, utility1, utility2)

if not(isfolder([pwd  '/data/energyL2' ]))                                           % create local directories for data storage
    mkdir([pwd  '/data/energyL2' ])
    mkdir([pwd  '/data/energyL2_t' ]);
    mkdir([pwd  '/data/spectrum' ]);
    mkdir([pwd  '/data/kappa' ]);
    mkdir([pwd  '/media' ]);
    mkdir([pwd  '/media/movies' ]);
    mkdir([pwd  '/media/energy' ]);
    mkdir([pwd  '/media/optimization' ]);
    mkdir([pwd  '/media/figures' ]);
    mkdir([pwd  '/media/figures/state' ]);
    mkdir([pwd  '/media/kappa' ]);
    addpath([pwd  '/data/energyL2' ])
    addpath([pwd  '/data/energyL2_t' ]);
    addpath([pwd  '/data/spectrum' ]);
    addpath([pwd  '/data/kappa' ]);
    addpath([pwd  '/media' ]);
    addpath([pwd  '/media/movies' ]);
    addpath([pwd  '/media/energy' ]);
    addpath([pwd  '/media/optimization' ]);
    addpath([pwd  '/media/figures' ]);
    addpath([pwd  '/media/figures/state' ]);
    addpath([pwd  '/media/kappa' ]);
end

parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
parfiglist = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];
optparfiglist = ['$\varphi = \tilde{\varphi}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];
originalIC = utility1;
tol = utility2(1);
if length(utility2) > 1
    if utility2(2) == 1
        optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) '_RCG' ];
    elseif utility2(2) == 0
        optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) '_RG' ];
    end
else
    optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) ];
end

% number of timesteps
time1 = ceil(T/dt); 
time2 = ceil(time1/save_each);
if T >= 0
    Ntime = max(time1,time2);
    Ntime_save = min(time1,time2);
    save_each = ceil(Ntime/Ntime_save);
    Ntime = ceil(Ntime/save_each);
end

%Ntime = size(u_n,2);
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

% 4-period physical space domain
x12_pts = 2*L_s1*linspace( 0 , 1 - 1/N , 2*N ); 
x22_pts = 2*L_s2*linspace( 0 , 1 - 1/N , 2*N ); 
[ x12x , x22x ] = meshgrid(x12_pts,x22_pts); % 2-dimensional grid

switch IC
    case {'optimized'}
        % compute L2 energy and Fourier mode evolution
        energyL2_og = NaN(Ntime,1);
        v_mean_og = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,Ntime);
        v_meancount_og = v_mean_og;
        Ntime_remaining = Ntime;
        for i = 1:Ntime

            if Ntime < Ntime_save_max && i == 1
                [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, T, N, K, L_s1, L_s2, Ntime, 0);
            else
                if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                    currentT = (i+Ntime_save_max-1)/Ntime*T;
                    [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, 0);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (mod(i,Ntime_save_max) == 1)
                    [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, T, N, K, L_s1, L_s2, Ntime_remaining, 0);
                end
            end
            
            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end
            u_i = reshape( u_og(:,imod) , [ N , N ] );
            energyL2_og(i,1) = (sum( u_og(:,imod) .* conj(u_og(:,imod)) )*(L1*L2)/N^2);
            v = fftshift(real(abs(fft2(u_i))));
            for j = 1:N
                for k = 1:N
                    index = round(sqrt((j-(N/2+1))^2+(k-(N/2+1))^2)) + 1;
                    v_mean_og(index,i) = v_mean_og(index,i) + v(j,k);
                    v_meancount_og(index,i) = v_meancount_og(index,i) + 1;
                end
            end
            for m = 1:size(v_meancount_og,1)
                v_mean_og(m,i) = v_mean_og(m,i)/v_meancount_og(m,i);
            end
                
            if i == 1
                u_IC_og = u_og(:,1);
            elseif i == Ntime
                u_TC_og = u_og(:,end);
            end

        end
        v_mean_og = v_mean_og(2:end,:);
        
        % L2 energy time derivative computation
        energyL2_t_og = NaN(Ntime-1,1);
        dt_save = T/(Ntime-1);
        for i = 2:length(energyL2_og)
            energyL2_t_og(i-1,1) = ( energyL2_og(i,1) - energyL2_og(i-1,1) ) / dt_save;
        end
end

switch solplot
    case {'norms','gif','diagnostics','initial','terminal','optdiag'}
        % compute L2 energy and Fourier mode evolution
        energyL2 = NaN(Ntime,1);
        v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,Ntime);
        v_meancount = v_mean;
        Ntime_remaining = Ntime;
        for i = 1:Ntime

            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                    currentT = (i+Ntime_save_max-1)/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (mod(i,Ntime_save_max) == 1)
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_remaining, utility1);
                end
            end
            
            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end
            u_i = reshape( u_n(:,imod) , [ N , N ] );
            energyL2(i,1) = (sum( u_n(:,imod) .* conj(u_n(:,imod)) )*(L1*L2)/N^2);
            v = fftshift(real(abs(fft2(u_i))));
            for j = 1:N
                for k = 1:N
                    index = round(sqrt((j-(N/2+1))^2+(k-(N/2+1))^2)) + 1;
                    v_mean(index,i) = v_mean(index,i) + v(j,k);
                    v_meancount(index,i) = v_meancount(index,i) + 1;
                end
            end
            for m = 1:size(v_meancount,1)
                v_mean(m,i) = v_mean(m,i)/v_meancount(m,i);
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
        
        energyL2data_file = [pwd '/data/energyL2/energyL2_' parameterlist '.dat'];
        energyL2deriv_file = [pwd '/data/energyL2_t/energyL2_' parameterlist '.dat'];
        spectrum_file = [pwd '/data/spectrum/spectrum_' parameterlist '.dat'];
        switch IC
            case {'optimized'}
                energyL2data_file = [pwd '/data/energyL2/energyL2_' optparameters '.dat'];
                energyL2deriv_file = [pwd '/data/energyL2_t/energyL2_' optparameters '.dat'];
                spectrum_file = [pwd '/data/spectrum/spectrum_' optparameters '.dat'];
        end
        writematrix(energyL2, energyL2data_file,'Delimiter','tab');
        writematrix(energyL2_t, energyL2deriv_file,'Delimiter','tab');
        writematrix(v_mean, spectrum_file,'Delimiter','tab');
end

switch solplot
    case 'gif'

<<<<<<< HEAD
        figure;
=======
        figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size
        

        set(gcf,'color','white')
        set(gca,'color','white')

        Ntime_remaining = Ntime;
        if utility2(end) > 0
            frames = utility2(end);
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

        for i = 1 : ceil(Ntime/frames) : Ntime+1
        
            if i == Ntime + 1
                i = Ntime;
            end
            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    currentT = frameovermax*Ntime_save_max/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (Ntime_remaining < Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_remaining, utility1);
                end
            end

            currentT = (i-1)/(Ntime-1)*T;

            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end

            subplot(2,2,1);
            % Draw surface plot
            u_i = reshape( u_n(:,imod) , [ N , N ] );
            surfc(x1pi,x2pi,u_i);
            xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
=======
            daspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
>>>>>>> a5f484f (Initial upload)
            view(3);
            drawnow

            subplot(2,2,2);
            % Draw surface plot
            u_i2x = [ u_i , u_i ; u_i, u_i];
            surfc(x12x,x22x,u_i2x);
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
=======
            %pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
>>>>>>> a5f484f (Initial upload)
            xline(L_s1,'--');
            yline(L_s2,'--');
            view(2);
            drawnow

            subplot(2,2,3);
            semilogy(timewindow,energyL2,'b')
            hold on
            xline(currentT,'-');
            hold off
            xlabel('Time $t$','Interpreter','latex');
            ylabel('$\| \phi(t;\varphi) \|^2_{L^2}$','Interpreter','latex');
            xlim([0 T])
            ylim([min(energyL2) max(energyL2)+1e-1])
            title("Evolution of $L^2$ energy",'Interpreter','latex')
            %legend('L^{2} energy','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean(:,i),"o--")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            title("Energy spectrum",'Interpreter','latex')
            xlim([1 size(v_mean,1)])
            ylim([ 1e-20 max(v_mean(1,:))+1e5 ])
            set(gca,'fontsize', 12) 
        
            title1 = 'Forward-time 2DKS solution';
            title2 = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            switch IC 
                case {'optimized'}
                    title2 = ['$\varphi = \tilde{\varphi}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            end
            sgtitle({title1, title2},'Interpreter','latex');

            if i == 1
<<<<<<< HEAD
                gif(filename)
=======
                gif(filename,'overwrite',true)
>>>>>>> a5f484f (Initial upload)
            else
                gif
            end
        
        end

        %{
        % manual conversion of GIF to AVI:  
        tic
        testdirs = dir([pwd '/media/movies/phys_sinL_N_48_T_2000_dt_0.01_Ls1_*.gif']);
        numberoftests = size(testdirs,1);
        for test = 1:numberoftests
            currenttest = testdirs(test).name;
            gif2avi(currenttest,[],'FrameRate',10)
            toc
        end
        %}

    case 'diagnostics'

        % Wavenumber evolution plot
        timewindow = linspace(0,T,Ntime);
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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

        %{
        % L2 energy time derivative plot
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        plot(timewindow(2:end),energyL2_t)
        set(gcf,'Position',[100 100 900 750])
        xlabel('Time $t$','Interpreter','latex'); 
        xlim([0 T])
        ylabel('$\frac{d}{dt} \| {\phi(t)} \|^2_{L^2}$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Evolution of $L^2$ energy time derivative','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        energyL2_t_file = [pwd '/media/energy/energyL2_deriv_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        exportgraphics(h,filename)
        %}

        switch IC 
            case {'optimized'}
                
                [diagnostics, linesearchJ] = load_2DKSsolution('optimization', 'optimized', dt, T, N, K, L_s1, L_s2, tol, originalIC);

                % L2 energy plot
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
                plotdata = diagnostics(:,7);
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD

                disp(['Optimization loop computation time: ' num2str(diagnostics(length(plotdata),5)) ])
=======
>>>>>>> a5f484f (Initial upload)
        end

    case 'initial'

<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])
        

        % inspect physical solution 
        u_T = reshape( u_IC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
=======
        %pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
>>>>>>> a5f484f (Initial upload)
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
        %{
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_initial'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters  '_initial'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])
        %}
        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_initial_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_initial_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

    case 'terminal'

<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])

        % inspect physical solution
        u_T = reshape( u_TC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
=======
        %pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
>>>>>>> a5f484f (Initial upload)
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
        %{
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_terminal'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])
        %}
        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_terminal_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

    case 'kappa'

        kappa = utility1;
        pertIC = utility2;
        kappaerror = abs( 1 - kappa );

<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        %epscheck = [.001,.0025,.005,.0075,.01,.025,.05,.075,.1,.25,.5,.75,1];
        epscheck = logspace(-15,-1,15);
        loglog(epscheck,kappaerror,'r-*')
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Kappa difference test','Interpreter','latex')
        subtitle(['$\varphi'' = \varphi_{' pertIC '}$, ' parfiglist],'Interpreter','latex','FontSize',14)
        filename = [pwd '/media/kappa/kappaerr_p' pertIC '_' parameterlist '.pdf'];
        exportgraphics(h,filename)

        kappa_file = [pwd '/data/kappa/kappa_p' pertIC '_' parameterlist '.dat'];
        writematrix(kappaerror, kappa_file,'Delimiter','tab');
    case 'optdiag'
        
        close all
        %% opt gif
<<<<<<< HEAD
        figure;
=======
        figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size
        

        set(gcf,'color','white')
        set(gca,'color','white')

        Ntime_remaining = Ntime;
        if utility2(end) > 0
            frames = utility2(end);
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

        for i = 1 : ceil(Ntime/frames) : Ntime+1
        
            if i == Ntime + 1
                i = Ntime;
            end
            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    currentT = frameovermax*Ntime_save_max/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (Ntime_remaining < Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_remaining, utility1);
                end
            end

            currentT = (i-1)/(Ntime-1)*T;

            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end

            subplot(2,2,1);
            % Draw surface plot
            u_i = reshape( u_n(:,imod) , [ N , N ] );
            surfc(x1pi,x2pi,u_i);
            xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
=======
            %pbaspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
>>>>>>> a5f484f (Initial upload)
            view(3);
            drawnow

            subplot(2,2,2);
            % Draw surface plot
            u_i2x = [ u_i , u_i ; u_i, u_i];
            surfc(x12x,x22x,u_i2x);
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
=======
            %pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
>>>>>>> a5f484f (Initial upload)
            xline(L_s1,'--');
            yline(L_s2,'--');
            view(2);
            drawnow

            subplot(2,2,3);
            semilogy(timewindow,energyL2,'b')
            hold on
            semilogy(timewindow,energyL2_og,'r--')
            xline(currentT,'-');
            hold off
            xlabel('Time $t$','Interpreter','latex');
            ylabel('$\| \phi(t;\varphi) \|^2_{L^2}$','Interpreter','latex');
            xlim([0 T])
            ylim([0.5*min([energyL2;energyL2_og]) 1.5*max([energyL2;energyL2_og]) ])
            title("Evolution of optimized $L^2$ energy",'Interpreter','latex')
            legend('$\tilde\varphi$','Interpreter','latex','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean(:,i),"o--")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            title("Energy spectrum",'Interpreter','latex')
            xlim([1 size(v_mean,1)])
            ylim([ 1e-20 1.5*max(max(v_mean)) ])
            set(gca,'fontsize', 12) 
        
            title1 = 'Forward-time 2DKS solution';
            title2 = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            switch IC 
                case {'optimized'}
                    title2 = ['$\varphi = \tilde{\varphi}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            end
            sgtitle({title1, title2},'Interpreter','latex');

            if i == 1
<<<<<<< HEAD
                gif(filename)
=======
                gif(filename,'overwrite',true)
>>>>>>> a5f484f (Initial upload)
            else
                gif
            end
        
        end

        close all
        %% diagnostics
        % Wavenumber evolution plot
        timewindow = linspace(0,T,Ntime);
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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

        switch IC 
            case {'optimized'}
                
                [diagnostics, linesearchJ] = load_2DKSsolution('optimization', 'optimized', dt, T, N, K, L_s1, L_s2, tol, originalIC);

                % L2 energy plot
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
                plotdata = diagnostics(:,7);
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD
                h = figure;
=======
                h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
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
<<<<<<< HEAD

                disp(['Optimization loop computation time: ' num2str(diagnostics(length(plotdata),5)) ])
        end

        %% opt initial
        h = figure;
=======
        end

        %% opt initial
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])
        

        % inspect physical solution 
        u_T = reshape( u_IC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
=======
        %pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
>>>>>>> a5f484f (Initial upload)
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
        %{
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_initial'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters  '_initial'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])
        %}
        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_initial_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_initial_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

        %% opt terminal
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])

        % inspect physical solution
        u_T = reshape( u_TC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
=======
        %pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
>>>>>>> a5f484f (Initial upload)
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
        %{
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_terminal'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])
        %}
        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_terminal_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

        %% finish opt plots
        IC = originalIC;
        parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
        parfiglist = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];

        %% orig initial
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])

        % inspect physical solution 
        u_T = reshape( u_IC_og , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
=======
        %pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
>>>>>>> a5f484f (Initial upload)
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
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_initial_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_initial_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

        %% orig terminal
<<<<<<< HEAD
        h = figure;
=======
        h = figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])

        % inspect physical solution
        u_T = reshape( u_TC_og , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
<<<<<<< HEAD
        pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
=======
        %pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
>>>>>>> a5f484f (Initial upload)
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
        % Save image
        filename = [pwd '/media/figures/state/phys_' parameterlist '_terminal_contour'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' optparameters '_terminal_contour'];
        end
        saveas(h,[filename '.fig'])
        exportgraphics(h,[filename '.pdf'])

        close all
        %% orig gif
<<<<<<< HEAD
        figure;
=======
        figure('Visible', 'off');
>>>>>>> a5f484f (Initial upload)
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size

        set(gcf,'color','white')
        set(gca,'color','white')

        Ntime_remaining = Ntime;
        if utility2(end) > 0
            frames = utility2(end);
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

        for i = 1 : ceil(Ntime/frames) : Ntime+1
        
            if i == Ntime + 1
                i = Ntime;
            end
            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    currentT = frameovermax*Ntime_save_max/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, K, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (Ntime_remaining < Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, K, L_s1, L_s2, Ntime_remaining, utility1);
                end
            end

            currentT = (i-1)/(Ntime-1)*T;

            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end

            subplot(2,2,1);
            % Draw surface plot
            u_i = reshape( u_n(:,imod) , [ N , N ] );
            surfc(x1pi,x2pi,u_i);
            xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
=======
            %pbaspect( [ abs(max(max(x1pi))), abs(max(max(x2pi))), abs(max(max(u_i))) ] );
>>>>>>> a5f484f (Initial upload)
            view(3);
            drawnow

            subplot(2,2,2);
            % Draw surface plot
            u_i2x = [ u_i , u_i ; u_i, u_i];
            surfc(x12x,x22x,u_i2x);
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
<<<<<<< HEAD
            pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
=======
            %pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
>>>>>>> a5f484f (Initial upload)
            xline(L_s1,'--');
            yline(L_s2,'--');
            view(2);
            drawnow

            subplot(2,2,3);
            semilogy(timewindow,energyL2_og,'b')
            hold on
            xline(currentT,'-');
            hold off
            xlabel('Time $t$','Interpreter','latex');
            ylabel('$\| \phi(t;\varphi) \|^2_{L^2}$','Interpreter','latex');
            xlim([0 T])
            ylim([0.5*min(energyL2_og) 1.5*max(energyL2_og) ])
            title("Evolution of $L^2$ energy",'Interpreter','latex')
            %legend('L^{2} energy','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean_og(:,i),"o--")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            title("Energy spectrum",'Interpreter','latex')
            xlim([1 size(v_mean_og,1)])
            ylim([ 1e-20 1.5*max(max(v_mean)) ])
            set(gca,'fontsize', 12) 
        
            title1 = 'Forward-time 2DKS solution';
            title2 = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            switch IC 
                case {'optimized'}
                    title2 = ['$\varphi = \tilde{\varphi}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(currentT,'%.2f') '$'];
            end
            sgtitle({title1, title2},'Interpreter','latex');

            if i == 1
<<<<<<< HEAD
                gif(filename)
=======
                gif(filename,'overwrite',true)
>>>>>>> a5f484f (Initial upload)
            else
                gif
            end
        
        end

end
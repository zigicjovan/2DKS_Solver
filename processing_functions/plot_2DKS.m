function plot_2DKS(save_each, solplot, IC, N, dt, T, L_s1, L_s2, Ntime_save_max, utility1, utility2)

% number of timesteps
time1 = ceil(T/dt); 
time2 = ceil(time1/save_each);
if T >= 1
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

% 4-period physical space domain
x12_pts = 2*L_s1*linspace( 0 , 1 - 1/N , 2*N ); 
x22_pts = 2*L_s2*linspace( 0 , 1 - 1/N , 2*N ); 
[ x12x , x22x ] = meshgrid(x12_pts,x22_pts); % 2-dimensional grid

switch IC
    case {'optimized'}
        originalIC = utility1;
        % compute L2 norm and Fourier mode evolution
        normL2_og = NaN(Ntime,1);
        v_mean_og = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,Ntime);
        v_meancount_og = v_mean_og;
        Ntime_remaining = Ntime;
        for i = 1:Ntime

            if Ntime < Ntime_save_max && i == 1
                [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, T, N, L_s1, L_s2, Ntime, 0);
            else
                if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                    currentT = (i+Ntime_save_max-1)/Ntime*T;
                    [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, currentT, N, L_s1, L_s2, Ntime_save_max, 0);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (mod(i,Ntime_save_max) == 1)
                    [u_og, ~] = load_2DKSsolution('forward', originalIC, dt, T, N, L_s1, L_s2, Ntime_remaining, 0);
                end
            end
            
            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end
            u_i = reshape( u_og(:,imod) , [ N , N ] );
            normL2_og(i,1) = sqrt(sum( u_og(:,imod) .* conj(u_og(:,imod)) )*(L1*L2)/N^2);
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
                %u_IC_og = u_og(:,1);
            elseif i == Ntime
                %u_TC_og = u_og(:,end);
            end

        end
        %v_mean_og = v_mean_og(2:end,:);
        
        % L2 norm time derivative computation
        normL2_t_og = NaN(Ntime-1,1);
        dt_save = T/(Ntime-1);
        for i = 2:length(normL2_og)
            normL2_t_og(i-1,1) = ( normL2_og(i,1) - normL2_og(i-1,1) ) / dt_save;
        end
end

switch solplot
    case {'norms','gif','diagnostics','initial','terminal'}
        % compute L2 norm and Fourier mode evolution
        normL2 = NaN(Ntime,1);
        v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,Ntime);
        v_meancount = v_mean;
        Ntime_remaining = Ntime;
        for i = 1:Ntime

            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (mod(i,Ntime_save_max) == 1)
                    currentT = (i+Ntime_save_max-1)/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (mod(i,Ntime_save_max) == 1)
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, L_s1, L_s2, Ntime_remaining, utility1);
                end
            end
            
            if mod(i,Ntime_save_max) ~= 0
                imod = mod(i,Ntime_save_max);
            else
                imod = Ntime_save_max;
            end
            u_i = reshape( u_n(:,imod) , [ N , N ] );
            normL2(i,1) = sqrt(sum( u_n(:,imod) .* conj(u_n(:,imod)) )*(L1*L2)/N^2);
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
        
        % L2 norm time derivative computation
        normL2_t = NaN(Ntime-1,1);
        dt_save = T/(Ntime-1);
        for i = 2:length(normL2)
            normL2_t(i-1,1) = ( normL2(i,1) - normL2(i-1,1) ) / dt_save;
        end
        
        mkdir([pwd  '/data/normL2' ]);
        mkdir([pwd  '/data/normL2_t' ]);
        mkdir([pwd  '/data/spectrum' ]);
        normL2data_file = [pwd '/data/normL2/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        normL2deriv_file = [pwd '/data/normL2_t/normL2_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        spectrum_file = [pwd '/data/spectrum/spectrum_' IC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        switch IC
            case {'optimized'}
                normL2data_file = [pwd '/data/normL2/normL2_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                normL2deriv_file = [pwd '/data/normL2_t/normL2_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
                spectrum_file = [pwd '/data/spectrum/spectrum_' IC '_' originalIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        end
        writematrix(normL2, normL2data_file,'Delimiter','tab');
        writematrix(normL2_t, normL2deriv_file,'Delimiter','tab');
        writematrix(v_mean, spectrum_file,'Delimiter','tab');
end

switch solplot
    case 'gif'

        figure;
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/media/movies' ]);

        set(gcf,'color','white')
        set(gca,'color','white')

        Ntime_remaining = Ntime;
        if utility2 > 0
            frames = utility2;
        else
            frames = Ntime;
        end
        frameovermax = 0;
        switch IC 
            case {'optimized'}
                frames = 100;
        end
        
        filename = [pwd '/media/movies/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_frames_' num2str(frames) '.gif'];

        switch IC
            case {'optimized'}
                filename = [pwd '/media/movies/phys_' IC '_' originalIC '_N_' num2str(N) ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_frames_' num2str(frames) '.gif'];
        end

        for i = 1 : ceil(Ntime/frames) : Ntime
        
            if Ntime < Ntime_save_max && i == 1
                [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, L_s1, L_s2, Ntime, utility1);
            else
                if (Ntime_remaining >= Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    currentT = frameovermax*Ntime_save_max/Ntime*T;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, currentT, N, L_s1, L_s2, Ntime_save_max, utility1);
                    Ntime_remaining = Ntime_remaining - Ntime_save_max;
                elseif (Ntime_remaining < Ntime_save_max) && (i > (frameovermax*Ntime_save_max))
                    frameovermax = frameovermax + 1;
                    [u_n, ~] = load_2DKSsolution('forward', IC, dt, T, N, L_s1, L_s2, Ntime_remaining, utility1);
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
            surfc(x1,x2,u_i);
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
            pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_i))) ] );
            view(3);
            drawnow

            subplot(2,2,2);
            % Draw surface plot
            u_i2x = [ u_i , u_i ; u_i, u_i];
            surfc(x12x,x22x,u_i2x);
            xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
            pbaspect( [ abs(max(max(x12x))), abs(max(max(x22x))), abs(max(max(u_i2x))) ] );
            xline(L_s1,'--');
            yline(L_s2,'--');
            view(2);
            drawnow

            subplot(2,2,3);
            semilogy(timewindow,normL2,'b')
            hold on
            xline(currentT,'-');
            hold off
            xlabel('Time $t$','Interpreter','latex');
            ylabel('$|| \phi(t;\varphi) ||$','Interpreter','latex');
            xlim([0 T])
            ylim([min(normL2) max(normL2)+1e-1])
            title("Evolution of $L^2$ norm",'Interpreter','latex')
            %legend('L^{2} norm','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean(:,i),".")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            title("Energy spectrum",'Interpreter','latex')
            xlim([1 size(v_mean,1)])
            ylim([1e-15 max(v_mean(1,:))])
            set(gca,'fontsize', 12) 
        
            title1 = 'Forward-time 2DKS solution';
            title2 = ['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(currentT,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'];
            switch IC 
                case {'optimized'}
                    title2 = ['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(currentT,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'];
            end
            sgtitle({title1, title2},'Interpreter','latex');

            if i == 1
                gif(filename)
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

        mkdir([pwd  '/media/energy' ]);

        % Wavenumber evolution plot
        timewindow = linspace(0,T,Ntime);
        h = figure;
        semilogy(timewindow,v_mean(1,:),'LineWidth',0.5,'Marker','.')
        hold on;
        for i = 2:size(v_mean,1)
            semilogy(timewindow,v_mean(i,:),'LineWidth',0.5,'Marker','.')
        end
        set(gcf,'Position',[100 100 900 750])
        xlabel('Time $t$','Interpreter','latex'); 
        xlim([0 T])
        ylim([1e-15 max(v_mean(1,:))+1e5 ])
        ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Evolution of Fourier spectrum','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        switch IC 
            case {'optimized'}
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        end
        %legend("Fourier band", 'Location','southeast','NumColumns',9,'Interpreter','latex')
        %frame = getframe(h);
        %im = frame2im(frame);
        filename = [pwd '/media/energy/wavenumberevol_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        filenamefig = [pwd '/media/energy/wavenumberevol_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/energy/wavenumberevol_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
                filenamefig = [pwd '/media/energy/wavenumberevol_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
        end
        %imwrite(im,filename,'png')
        exportgraphics(h,filename)
        saveas(h,filenamefig)

        % L2 norm plot
        h = figure;
        semilogy(timewindow,normL2,'LineWidth',0.5,'Marker','.')
        set(gcf,'Position',[100 100 900 750])
        xlabel('Time $t$','Interpreter','latex'); 
        xlim([0 T])
        ylabel('$||{\phi(t;\varphi)}||_{L^2}$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Evolution of $L^2$ norm','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        switch IC 
            case {'optimized'}
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        end
        filename = [pwd '/media/energy/normL2_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        filenamefig = [pwd '/media/energy/normL2_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/energy/normL2_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/energy/normL2_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
        end
        exportgraphics(h,filename)
        saveas(h,filenamefig)

        %{
        % L2 norm time derivative plot
        h = figure;
        plot(timewindow(2:end),normL2_t)
        set(gcf,'Position',[100 100 900 750])
        xlabel('Time $t$','Interpreter','latex'); 
        xlim([0 T])
        ylabel('$\frac{d}{dt} ||{\phi(t)}||_{L^2}$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Evolution of $L^2$ norm time derivative','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        normL2_t_file = [pwd '/media/energy/normL2_deriv_' IC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        exportgraphics(h,filename)
        %}

        switch IC 
            case {'optimized'}
                mkdir([pwd  '/media/optimization' ]);
                [diagnostics, linesearchJ] = load_2DKSsolution('optimization', 'optimized', dt, T, N, L_s1, L_s2, 0, originalIC);

                % L2 norm plot
                h = figure;
                semilogy(timewindow,normL2_og,'LineWidth',0.5,'Marker','.')
                hold on
                semilogy(timewindow,normL2,'LineWidth',0.5,'Marker','.')
                hold off
                set(gcf,'Position',[100 100 900 750])
                xlabel('Time $t$','Interpreter','latex'); 
                xlim([0 T])
                ylabel('$||{\phi(t;\varphi)}||_{L^2}$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Evolution of optimized $L^2$ norm','Interpreter','latex')
                subtitle(['$L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                legend('$\phi(t,\varphi^{(0)})$','$\phi(t,\widetilde{\varphi})$','Interpreter','latex','Location','northwest')
                filename = [pwd '/media/optimization/normL2comp_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/normL2comp_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % J history plot
                plotdata = diagnostics(:,1);
                h = figure;
                plot(plotdata,'r-*')
                set(gcf,'Position',[100 100 900 750])
                xlabel('Iteration number $n$','Interpreter','latex'); 
                xlim([1 length(plotdata)])
                ylabel('$\mathcal{J}_T(\varphi^{(n)})$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Evolution of objective functional $\mathcal{J}_T(\varphi^{(n)})=||{\phi^{(n)}(T)}||^2_{L^2}$','Interpreter','latex')
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/Jhistory_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/Jhistory_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % J change plot
                plotdata = diagnostics(:,2);
                h = figure;
                semilogy(plotdata,'r-*')
                set(gcf,'Position',[100 100 900 750])
                xlabel('Iteration number $n$','Interpreter','latex'); 
                xlim([1 length(plotdata)])
                ylabel('$\frac{\mathcal{J}_T(\varphi^{(n+1)}) - \mathcal{J}_T(\varphi^{(n)})}{\mathcal{J}_T(\varphi^{(n)})}$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Relative change in objective functional $\mathcal{J}_T(\varphi^{(n)})=||{\phi^{(n)}(T)}||^2_{L^2}$','Interpreter','latex')
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/Jchange_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/Jchange_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % stepsize plot
                plotdata = diagnostics(:,3);
                h = figure;
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
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/stepsize_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/stepsize_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % manifold plot
                plotdata = diagnostics(:,4);
                h = figure;
                semilogy(plotdata,'r-*')
                set(gcf,'Position',[100 100 900 750])
                xlabel('Iteration number $n$','Interpreter','latex'); 
                xlim([1 length(plotdata)])
                ylabel('$||\varphi^{(n)}||^2_{L^2}$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Evolution of initial condition magnitude $K=||\varphi^{(n)}||^2_{L^2}$','Interpreter','latex')
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/initialK_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/initialK_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % gradient magnitude plot
                plotdata = diagnostics(:,6);
                h = figure;
                semilogy(plotdata,'r-*')
                set(gcf,'Position',[100 100 900 750])
                xlabel('Iteration number $n$','Interpreter','latex'); 
                xlim([1 length(plotdata)])
                ylabel('$||\nabla \mathcal{J}_T(\varphi^{(n)})||^2_{L^2}$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Evolution of objective gradient magnitude $||\nabla \mathcal{J}_T(\varphi^{(n)})||^2_{L^2}$','Interpreter','latex')
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/Jgradient_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/Jgradient_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % momentum magnitude plot
                plotdata = diagnostics(:,7);
                h = figure;
                semilogy(plotdata,'r-*')
                set(gcf,'Position',[100 100 900 750])
                xlabel('Iteration number $n$','Interpreter','latex'); 
                xlim([1 length(plotdata)])
                ylabel('$||\beta^{(n)}||^2_{L^2}$','Interpreter','latex');
                fontsize(12,"points")
                set(gca,'fontsize', 16) 
                set(gcf,'color','white')
                set(gca,'color','white')    
                title('Evolution of momentum magnitude $||\beta^{(n)}||^2_{L^2}$','Interpreter','latex')
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/momentumsize_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/momentumsize_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                % line search plot
                h = figure;
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
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
                filename = [pwd '/media/optimization/linesearch_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
                filenamefig = [pwd '/media/optimization/linesearch_' IC '_' originalIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.fig'];
                exportgraphics(h,filename)
                saveas(h,filenamefig)

                close all

                disp(['Optimization loop computation time: ' num2str(diagnostics(length(plotdata),5)) ])
        end

    case 'initial'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])
        mkdir([pwd  '/media/figures/state' ]);

        % inspect physical solution 
        u_T = reshape( u_IC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
        pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white') 
        colormap(redblue)
        view(3);
        title('Initial state','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        switch IC 
            case {'optimized'}
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        end
        % Save image
        filename = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial.pdf'];
        filenamefig = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial.pdf'];
                filenamefig = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial.fig'];
        end
        exportgraphics(h,filename)
        saveas(h,filenamefig)

        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial_contour.pdf'];
        filenamefig = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial_contour.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial_contour.pdf'];
                filenamefig = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_initial_contour.fig'];
        end
        exportgraphics(h,filename)
        saveas(h,filenamefig)

    case 'terminal'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        set(gcf,'Position',[100 100 900 750])
        mkdir([pwd  '/media/figures/state' ]);

        % inspect physical solution
        u_T = reshape( u_TC , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('$\frac{x_1}{2\pi}$','Interpreter','latex'); ylabel('$\frac{x_2}{2\pi}$','Interpreter','latex'); zlabel('u(x_1,x_2)');
        shading interp
        pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
        view(3);
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white') 
        colormap(redblue)
        title('Terminal state','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{' IC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        switch IC 
            case {'optimized'}
                subtitle(['$\varphi = \widetilde{\varphi}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.1f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        end
        % Save image
        filename = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal.pdf'];
        filenamefig = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal.pdf'];
                filenamefig = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal.fig'];
        end
        exportgraphics(h,filename)
        saveas(h,filenamefig)

        % contour plot
        view(2);
        % Save image
        filename = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal_contour.pdf'];
        filenamefig = [pwd '/media/figures/state/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal_contour.fig'];
        switch IC
            case {'optimized'}
                filename = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal_contour.pdf'];
                filenamefig = [pwd '/media/figures/state/phys_' IC '_' originalIC '_N_' num2str(N) ...
                    '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '_terminal_contour.fig'];
        end
        exportgraphics(h,filename)
        saveas(h,filenamefig)

    case 'kappa'

        mkdir([pwd  '/media/kappa' ]);

        kappa = utility1;
        pertIC = utility2;
        kappaerror = abs( 1 - kappa );

        h = figure;
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
        subtitle(['$\varphi = \varphi_{' IC '}, \varphi'' = \varphi_{' pertIC '}, L_1 = 2\pi(' num2str(L_s1,'%.3f') '), L_2 = 2\pi(' num2str(L_s2,'%.3f') '), T = ' num2str(T,'%.0f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        filename = [pwd '/media/kappa/kappaerr_' IC '_p' pertIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        exportgraphics(h,filename)

        mkdir([pwd  '/data/kappa' ]);
        kappa_file = [pwd '/data/kappa/kappa_' IC '_p' pertIC '_N_' num2str(N) '' ...
                '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        writematrix(kappaerror, kappa_file,'Delimiter','tab');

end
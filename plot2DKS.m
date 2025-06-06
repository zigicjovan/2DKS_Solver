function plot2DKS(u_n, solplot, IC, N, dt, T, L_s1, L_s2, utility)

Ntime = size(u_n,2);
timewindow = linspace(0,T,Ntime);

% length-scale parameters
L1 = 2*pi*L_s1;
L2 = 2*pi*L_s2;

% unit physical space domain
x1_pts = L1*linspace( 0 , 1 - 1/N , N ); 
x2_pts = L2*linspace( 0 , 1 - 1/N , N ); 
[ x1 , x2 ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

% repeated physical space domain
x12_pts = 2*L1*linspace( 0 , 1 - 1/N , 2*N ); 
x22_pts = 2*L2*linspace( 0 , 1 - 1/N , 2*N ); 
[ x12x , x22x ] = meshgrid(x12_pts,x22_pts); % 2-dimensional grid

% fourier space domain for nonlinear term
k1_0_pts = [ 0 : N/2-1 , 0 , -N/2+1 : -1]; 
k2_0_pts = [ 0 : N/2-1 , 0 , -N/2+1 : -1]; 
[ k1_0 , k2_0 ] = meshgrid(k1_0_pts,k2_0_pts); % 2-dimensional grid

% compute L2 norm and Fourier mode evolution
normL2 = NaN(Ntime,1);
v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,Ntime);
v_meancount = v_mean;
for i = 1:Ntime
    u_i = reshape( u_n(:,i) , [ N , N ] );
    normL2(i,1) = sqrt(sum( u_n(:,i) .* conj(u_n(:,i)) )*(L1*L2)/N^2);
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
end
v_mean = v_mean(2:end,:);

% L2 norm time derivative computation
normL2_t = NaN(Ntime-1,1);
dt_save = T/(Ntime-1);
for i = 2:length(normL2)
    normL2_t(i-1,1) = ( normL2(i,1) - normL2(i-1,1) ) / dt_save;
end

switch solplot
    case 'gif'

        figure;
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/movies' ]);
        filename = [pwd '/data/media/movies/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '.gif'];
        set(gcf,'color','white')
        set(gca,'color','white')
        
        for i = 1 : ceil(Ntime/Ntime) : Ntime
        
            currentT = (i-1)/(Ntime-1)*T;

            subplot(2,2,1);
            % Draw surface plot
            u_i = reshape( u_n(:,i) , [ N , N ] );
            surfc(x1,x2,u_i);
            xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
            pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_i))) ] );
            view([37.5,30]);
            drawnow

            subplot(2,2,2);
            % Draw surface plot
            u_i2x = [ u_i , u_i ; u_i, u_i];
            surfc(x12x,x22x,u_i2x);
            xlabel('x_1'); ylabel('x_2'); %zlabel('Solution')
            shading(gca,'interp')
            colormap(redblue)
            pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_i)) ] );
            xline(L1,'--');
            yline(L2,'--');
            view(2);
            drawnow

            subplot(2,2,3);
            semilogy(timewindow,normL2,'b')
            hold on
            xline(currentT,'-');
            hold off
            xlabel('$t$','Interpreter','latex');
            ylabel('$|| \phi(t) ||$','Interpreter','latex');
            xlim([0 T])
            ylim([min(normL2) max(normL2)+1e-1])
            title("Evolution of $L^2$ norm",'Interpreter','latex')
            legend('L^{2} norm','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean(:,i),".")
            xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex'); 
            ylabel('$\frac{1}{j}\sum_{j} |{\widehat\phi_k}|$','Interpreter','latex');
            title("Energy spectrum")
            xlim([1 size(v_mean,1)])
            ylim([1e-15 max(v_mean(1,:))])
            set(gca,'fontsize', 12) 
        
            sgtitle(['2DKS, x_1 = 2\pi(' num2str(L_s1) '), x_2 = 2\pi(' num2str(L_s2) '), T = ' num2str(currentT,'%.2f') ', {\Delta}t = ' num2str(dt) ', N = ' num2str(N) ''])
        
            if i == 1
                gif(filename)
            else
                gif
            end
        
        end

    case 'diagnostics'

         mkdir([pwd  '/data/energy' ]);

        % Wavenumber evolution plot
        timewindow = linspace(0,T,Ntime);
        h = figure;
        semilogy(timewindow,v_mean(1,:))
        hold on;
        for i = 2:size(v_mean,1)
            semilogy(timewindow,v_mean(i,:))
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
        title("Evolution of Fourier spectrum")
        legend("Fourier mode", 'Location','southeast','NumColumns',9,'Interpreter','latex')
        frame = getframe(h);
        im = frame2im(frame);
        wavenumberevol_file = [pwd '/data/energy/wavenumberevol_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,wavenumberevol_file,'png');

        % L2 norm plot
        h = figure;
        semilogy(timewindow,normL2)
        set(gcf,'Position',[100 100 900 750])
        xlabel('Time $t$','Interpreter','latex'); 
        xlim([0 T])
        ylabel('$||{\phi(t)}||_{L^2}$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Evolution of L2 norm")
        frame = getframe(h);
        im = frame2im(frame);
        normL2_file = [pwd '/data/energy/normL2_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,normL2_file,'png');

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
        title("Evolution of L2 norm time derivative")
        frame = getframe(h);
        im = frame2im(frame);
        normL2_t_file = [pwd '/data/energy/normL2_deriv_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,normL2_t_file,'png');

    case 'initial'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/figures' ]);

        % inspect physical and fourier spectrum
        v_T = reshape( fft2(u_n(:,1)) , [ N , N ] ); 
        u_T = reshape( u_n(:,1) , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('x_1'); ylabel('x_2'); zlabel('u(x_1,x_2)');
        shading interp
        pbaspect( [ max(max(x1)), max(max(x2)), max(max(u_T)) ] );
        set(gcf,'color','white')
        set(gca,'color','white')
        colormap(redblue)
        view(3);
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_initial.png'];
        imwrite(im,filename,'png');

        % contour plot
        view(2);
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_initial_contour.png'];
        imwrite(im,filename,'png');

        % fourier plot
        surfc(k1_0,k2_0,abs(v_T)); set(gca,'xscale','log','yscale','log','zscale','log');
        xlabel('k_1'); ylabel('k_2'); zlabel('|v(k_1,k_2)|');
        set(gcf,'color','white')
        set(gca,'color','white')
        colormap(redblue)
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/four_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_initial.png'];
        imwrite(im,filename,'png');

    case 'terminal'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/figures' ]);

        % inspect physical and fourier spectrum
        v_T = reshape( fft2(u_n(:,Ntime)) , [ N , N ] ); 
        u_T = reshape( u_n(:,Ntime) , [ N , N ] );

        % surface plot
        surfc(x1,x2,u_T); 
        xlabel('x_1'); ylabel('x_2'); zlabel('u(x_1,x_2)');
        shading interp
        pbaspect( [ abs(max(max(x1))), abs(max(max(x2))), abs(max(max(u_T)))] );
        view(3);
        set(gcf,'color','white')
        set(gca,'color','white')
        colormap(redblue)
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_terminal.png'];
        imwrite(im,filename,'png');

        % contour plot
        view(2);
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_terminal_contour.png'];
        imwrite(im,filename,'png');

        % fourier plot
        surfc(k1_0,k2_0,real(abs(v_T))); set(gca,'xscale','log','yscale','log','zscale','log');
        xlabel('k_1'); ylabel('k_2'); zlabel('|v(k_1,k_2)|');
        set(gcf,'color','white')
        set(gca,'color','white')
        colormap(redblue)
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/four_' IC '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_terminal.png'];
        imwrite(im,filename,'png');

    case 'kappa'

        mkdir([pwd  '/data/kappa' ]);

        kappa = utility;
        kappaerror = abs( 1 - kappa );

        h = figure;
        semilogx(logspace(-15,-1,15),kappa)
        yline(1,'--')
        %set(gcf,'Position',[100 100 900 750])
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylim([0.97 1.03])
        xlim([1e-15 1e-1])
        ylabel('$\kappa(\varepsilon)$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa test for $\varphi$", 'Interpreter','latex')
        frame = getframe(h);
        im = frame2im(frame);
        kappa1_file = [pwd '/data/kappa/kappa1_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappa1_file,'png');

        h = figure;
        loglog(logspace(-15,-1,15),kappaerror)
        %set(gcf,'Position',[100 100 900 750])
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa difference test for $\varphi_1$", 'Interpreter','latex')
        frame = getframe(h);
        im = frame2im(frame);
        kappaerr_file = [pwd '/data/kappa/kappaerr_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappaerr_file,'png');

        %{
        kappasinL = [ kappalist_48_1e2_sinL_p1 kappalist_48_5e3_sinL_p1 kappalist_48_1e3_sinL_p1 ...
        kappalist_64_1e2_sinL_p1 kappalist_64_5e3_sinL_p1 kappalist_64_1e3_sinL_p1 ] ;
        kappasin = [ kappalist_48_1e2_sin_p1 kappalist_48_5e3_sin_p1 kappalist_48_1e3_sin_p1 ...
        kappalist_64_1e2_sin_p1 kappalist_64_5e3_sin_p1 kappalist_64_1e3_sin_p1 ] ;

        kappaerrorsinL= abs( 1 - kappasinL );
        kappaerrorsin = abs( 1 - kappasin );

        %%% first test %%%
        h = figure;
        semilogx(logspace(-15,-1,15),kappasinL(:,1),'r--+')
        hold on
        semilogx(logspace(-15,-1,15),kappasinL(:,4),'r--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,7),'r--o')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+1),'r:+')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+4),'r:*')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+7),'r:o')
        semilogx(logspace(-15,-1,15),kappasinL(:,2),'g--+')
        semilogx(logspace(-15,-1,15),kappasinL(:,5),'g--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,8),'g--o')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+2),'g:+')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+5),'g:*')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+8),'g:o')
        semilogx(logspace(-15,-1,15),kappasinL(:,3),'b--+')
        semilogx(logspace(-15,-1,15),kappasinL(:,6),'b--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,9),'b--o')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+3),'b:+')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+6),'b:*')
        %semilogx(logspace(-15,-1,15),kappasinL(:,9+9),'b:o')

        yline(1,'--')
        %set(gcf,'Position',[100 100 900 750])
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylim([0.97 1.03])
        xlim([1e-15 1e-1])
        ylabel('$\kappa(\varepsilon)$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa test for $\varphi_1$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa1_file = [pwd '/data/kappa/kappasinL_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappa1_file,'png');

        h = figure;
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,1),'r--+')
        hold on
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,4),'r--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,7),'r--o')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+1),'r:+')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+4),'r:*')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+7),'r:o')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,2),'g--+')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,5),'g--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,8),'g--o')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+2),'g:+')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+5),'g:*')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+8),'g:o')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,3),'b--+')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,6),'b--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,9),'b--o')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+3),'b:+')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+6),'b:*')
        %loglog(logspace(-15,-1,15),kappaerrorsinL(:,9+9),'b:o')

        %set(gcf,'Position',[100 100 900 750])
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa difference test for $\varphi_1$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsinL_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappa2_file,'png');

        %%% second test %%%
        h = figure;
        semilogx(logspace(-15,-1,15),kappasin(:,1),'r--+')
        hold on
        semilogx(logspace(-15,-1,15),kappasin(:,4),'r--*')
        semilogx(logspace(-15,-1,15),kappasin(:,7),'r--o')
        semilogx(logspace(-15,-1,15),kappasin(:,2),'g--+')
        semilogx(logspace(-15,-1,15),kappasin(:,5),'g--*')
        semilogx(logspace(-15,-1,15),kappasin(:,8),'g--o')
        semilogx(logspace(-15,-1,15),kappasin(:,3),'b--+')
        semilogx(logspace(-15,-1,15),kappasin(:,6),'b--*')
        semilogx(logspace(-15,-1,15),kappasin(:,9),'b--o')

        yline(1,'--')
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        xlim([1e-15 1e-1])
        ylim([0.97 1.03])
        ylabel('$\kappa(\varepsilon)$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa test for $\varphi_2$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa1_file = [pwd '/data/kappa/kappasin_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappa1_file,'png');

        h = figure;
        loglog(logspace(-15,-1,15),kappaerrorsin(:,1),'r--+')
        hold on
        loglog(logspace(-15,-1,15),kappaerrorsin(:,4),'r--*')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,7),'r--o')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,2),'g--+')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,5),'g--*')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,8),'g--o')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,3),'b--+')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,6),'b--*')
        loglog(logspace(-15,-1,15),kappaerrorsin(:,9),'b--o')

        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa difference test for $\varphi_2$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsin_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
        imwrite(im,kappa2_file,'png');
        %}

end
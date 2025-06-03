function plot2DKS(v_n , u_n, solplot, method, N, dt, T, L_s1, L_s2, normL2, v_mean)

Ntime = size(u_n,2);
timewindow = linspace(0,T,Ntime);

% length-scale parameters
L_x1 = (1/L_s1);
L_x2 = (1/L_s2);
L1 = 2*pi/L_x1;
L2 = 2*pi/L_x2;

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

switch solplot
    case 'gif'

        figure;
        set(gcf,'Position',[100 100 900 750])
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/movies' ]);
        filename = [pwd '/data/media/movies/phys_' method '_N_' num2str(N) ...
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
            semilogy(timewindow,normL2,'b','LineWidth',1)
            hold on
            xline(currentT,'-','LineWidth',1);
            hold off
            xlabel('$t$','Interpreter','latex');
            ylabel('$|| \phi(t) ||$','Interpreter','latex');
            xlim([0 T])
            ylim([min(normL2) max(normL2)+1e-1])
            title("Evolution of $L^2$ norm",'Interpreter','latex')
            legend('L^{2} norm','Location','southeast')
            set(gca,'fontsize', 12) 
        
            subplot(2,2,4);
            semilogy(v_mean(:,i),".",'LineWidth',1)
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

    case 'timestrip'
        
        % Plot strip of solution over time
        line = 1; 
        strips = 1;
        N_xs = N*strips;
        u_l = NaN( N_xs , Ntime );
        x1s = size( x1 , size(x1,2) * strips );
        x2s = size( x2 , size(x2,2) * strips );
        x1s( x1 , 1:size(x1,2) ) = x1;
        x2s( x2 , 1:size(x2,2) ) = x2;
        for i = 2:strips
            col_start_x1 = (i-1) * size(x1,2) + 1;
            col_end_x1 = (i) * size(x1,2);
            col_start_x2 = (i-1) * size(x2,2) + 1;
            col_end_x2 = (i) * size(x2,2);
            x1s( x1 , col_start_x1:col_end_x1 ) = [ x1s ; x1 ];
            x2s( x2 , col_start_x2:col_end_x2 ) = [ x2s ; x2 ];
        end
        for i = 1:Ntime
            u_l1 = reshape( u_n( line*N + 1 : (strips+1)*N , i ) , [ N_xs , 1 ] );
            u_l( : , i ) = u_l1;
        end
        figure(2);
        surf( x1s , x2s , u_l )
        shading interp

    case 'initial'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/figures' ]);

        % inspect physical and fourier spectrum
        v_T = reshape( v_n(:,1) , [ N , N ] ); 
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
        filename = [pwd '/data/media/figures/phys_' method '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_initial.png'];
        imwrite(im,filename,'png');

        % contour plot
        view(2);
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' method '_N_' num2str(N) ...
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
        filename = [pwd '/data/media/figures/four_' method '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_initial.png'];
        imwrite(im,filename,'png');

    case 'terminal'

        h = figure;
        axis tight manual % this ensures that getframe() returns a consistent size
        mkdir([pwd  '/data/media/figures' ]);

        % inspect physical and fourier spectrum
        v_T = reshape( v_n(:,Ntime) , [ N , N ] ); 
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
        filename = [pwd '/data/media/figures/phys_' method '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_terminal.png'];
        imwrite(im,filename,'png');

        % contour plot
        view(2);
        % Save image
        frame = getframe(h);
        im = frame2im(frame);
        filename = [pwd '/data/media/figures/phys_' method '_N_' num2str(N) ...
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
        filename = [pwd '/data/media/figures/four_' method '_N_' num2str(N) ...
            '_T_' num2str(T) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1) '_Ls2_' num2str(L_s2) '_terminal.png'];
        imwrite(im,filename,'png');

end
function [measure1,measure2,measure3] = plot_measures(testcase, timestep, control, gridsize, T, L_s1, utility1, utility2)

switch testcase
    case 'kappa'

        mkdir([pwd  '/media/kappa' ]);
        measure1 = 0;
        measure2 = 0;
        measure3 = 0;

        start = utility1 - utility2;
        dt = timestep;
        IC = control;
        pertIC = control;
        N = gridsize;
        lastT = T(utility2);
        L_s2 = L_s1;

        kappalist_file = [pwd '/data/kappa/kappalist_' IC '_p' pertIC '_N_' num2str(N) '' ...
            '_T_' num2str(lastT) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];

        rieszlist_file = [pwd '/data/kappa/rieszlist_' IC '_p' pertIC '_N_' num2str(N) '' ...
             '_T_' num2str(lastT) '_dt_' num2str(dt) '_Ls1_' num2str(L_s1,'%.3f') '_Ls2_' num2str(L_s2,'%.3f') '.dat'];
        
        kappalist = readmatrix(kappalist_file);
        rieszlist = readmatrix(rieszlist_file);

        kappalist = abs(1 - kappalist);

        h = figure;
        semilogy(T,rieszlist(start+1:start+utility2,1),'r-*')
        title('Kappa test denominator','Interpreter','latex')
        subtitle(['$\ell = ' num2str(L_s1,'%.3f') ', \varphi_{' IC '}, \Delta t = ' num2str(dt,'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        xlabel('Time window $T$','Interpreter','latex'); 
        ylabel('$\langle \nabla\mathcal{J}_T (\varphi), \varphi''  \rangle_{L^2}$','Interpreter','latex');
        xlim([ T(1) T(utility2) ])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')        
        legend(['$\ell = ' num2str(L_s1,'%.3f') ', \varphi_{' IC '}$'],...
            'Interpreter','latex','Location','northwest')
        hold off
        filename = [pwd '/media/kappa/rieszlist_' IC '_p' pertIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(lastT) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        exportgraphics(h,filename)

        lgd = cell(utility2,1);
        h = figure;
        loglog(logspace(-15,-1,15),kappalist(:,start+1),'r--x')
        lgd{1} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(1),'%.0f') '$'];
        hold on
        loglog(logspace(-15,-1,15),kappalist(:,start+2),'r--d')
        lgd{2} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(2),'%.0f') '$'];
        loglog(logspace(-15,-1,15),kappalist(:,start+3),'r--o')
        lgd{3} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(3),'%.0f') '$'];
        loglog(logspace(-15,-1,15),kappalist(:,start+4),'g--x')
        lgd{4} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(4),'%.0f') '$'];
        loglog(logspace(-15,-1,15),kappalist(:,start+5),'g--d')
        lgd{5} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(5),'%.0f') '$'];
        if utility2 == 10
            loglog(logspace(-15,-1,15),kappalist(:,start+6),'g--o')
            lgd{6} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(6),'%.0f') '$'];
            loglog(logspace(-15,-1,15),kappalist(:,start+7),'b--x')
            lgd{7} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(7),'%.0f') '$'];
            loglog(logspace(-15,-1,15),kappalist(:,start+8),'b--d')
            lgd{8} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(8),'%.0f') '$'];
            loglog(logspace(-15,-1,15),kappalist(:,start+9),'b--o')
            lgd{9} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(9),'%.0f') '$'];
            loglog(logspace(-15,-1,15),kappalist(:,start+10),'m--x')
            lgd{10} = ['$\ell = ' num2str(L_s1,'%.3f') ', T = ' num2str(T(10),'%.0f') '$'];
        end
        title('Kappa difference test','Interpreter','latex')
        subtitle(['$\ell = ' num2str(L_s1,'%.3f') ', \varphi_{' IC '}, \Delta t = ' num2str(dt,'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex','FontSize',14)
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')        
        legend(lgd,'Interpreter','latex','Location','north','NumColumns',3)
        hold off
        filename = [pwd '/media/kappa/kappalist_' IC '_p' pertIC '_N' num2str(N) '_dt' num2str(dt) '_T' num2str(lastT) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.pdf'];
        exportgraphics(h,filename)

        close all

    case 'energygrowth'

        IC = control;
        L_scale = timestep;
        L_target = L_s1;
        l2norms_mode  = gridsize;
        l2norms_avg = utility1;

        figure
        plot(L_scale,l2norms_mode(:,1),'r-x')
        hold on
        try
            plot(L_scale,l2norms_mode(:,2),'b-x')
            plot(L_scale,l2norms_mode(:,3),'g-x')
            plot(L_scale,l2norms_mode(:,4),'k-x')
        catch
        end
        plot(L_scale,l2norms_avg(:,1),'r-o')
        try
            plot(L_scale,l2norms_avg(:,2),'b-o')
            plot(L_scale,l2norms_avg(:,3),'g-o')
            plot(L_scale,l2norms_avg(:,4),'k-o')
        catch
        end
        for i = 1:length(L_target)
            xl = xline(L_target(i),'--',{num2str(L_target(i),'%.3f')});
            xl.LabelVerticalAlignment = 'top';
            xl.LabelHorizontalAlignment = 'left';
        end
        title(['Modal energy ($t\in [0,' num2str(T) ']$) and mean energy ($t\in [' num2str(T/2) ',' num2str(T) ']$)'],'Interpreter','latex')
        xlabel('Domain factor $\ell$','Interpreter','latex')
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white') 
        legendcell = cell(1,2*length(IC));
        for i = 1:length(IC)
            legendcell(i) = {['$\check{\phi}^{L^2}_{[0,' num2str(T) ']}(\varphi_{' strjoin(IC(i),'') '})$']};
            legendcell(length(IC) + i) = {['$\bar{\phi}^{L^2}_{[' num2str(T/2) ',' num2str(T) ']}(\varphi_{' strjoin(IC(i),'') '})$']};
        end
        legend(legendcell,'Location','northwest','Interpreter','latex')
        legend('boxoff')
        hold off

        figure
        plot(L_scale,l2norms_mode(:,end),'r-*')
        hold on
        plot(L_scale,l2norms_avg(:,end),'b-o')
        for i = 1:length(L_target)
            xl = xline(L_target(i),'--',{num2str(L_target(i),'%.3f')});
            xl.LabelVerticalAlignment = 'top';
            xl.LabelHorizontalAlignment = 'left';
        end
        title(['Difference in values for ${IC} = \{ ' strjoin(IC(1),'') ',' strjoin(IC(2),'') ' \}$ within $T = ' num2str(T) '$'],'Interpreter','latex')
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')  
        xlabel('Domain factor $\ell$','Interpreter','latex')
        legendcell = cell(1,2);
        legendcell(1) = {['$| \check{\phi}^{L^2}_{[0,' num2str(T) ']}(\varphi_{' strjoin(IC(1),'') '}) - \check{\phi}^{L^2}_{[0,' num2str(T) ']}(\varphi_{' strjoin(IC(2),'') '}) |$' ]};
        legendcell(2) = {['$| \bar{\phi}^{L^2}_{[' num2str(T/2) ',' num2str(T) ']}(\varphi_{' strjoin(IC(1),'') '}) - \bar{\phi}^{L^2}_{[' num2str(T/2) ',' num2str(T) ']}(\varphi_{' strjoin(IC(2),'') '}) |$' ]};
        legend(legendcell,'Location','northwest','Interpreter','latex')

        measure1 = 0;
        measure2 = 0;
        measure3 = 0;
    case 'temporal'
        L_s2 = utility1;
        %%% (1) choose smallest step %%%
        dt = timestep(end); % size of smallest timestep
    
        %%% (2) load reference solution %%%
    
        [u_n, ~, time_n] = load_2DKSsolution('time_evolution', control, dt, T, gridsize, L_s1, L_s2);
        u_Nsq = reshape( u_n(:,end) , gridsize, gridsize);
        measure1 = NaN( length(timestep) , 2 );
        measure2 = measure1;
        measure3 = NaN( length(timestep) , 2 );
    
        %%% (3) error analysis %%%
        measure1( end , 1) = dt;
        measure2( end , 1) = dt;
        measure1( end , 2) = 0;
        measure2( end , 2) = 0;
    
        %%% (4) computational time %%%
        measure3( end , 1) = dt;
        measure3( end , 2) = time_n;
    
        j = 1; % iteration counter for saved measures
    
        for i = 1 : length(timestep)-1 % plot linf, l2 errors againt dt_small
    
            %%% (1) choose larger step %%%
            dt = timestep(i); % size of comparison timestep
    
            %%% (2) load solution %%%
            [u_n, ~, time_n] = load_2DKSsolution('time_evolution', control, dt, T, gridsize, L_s1, L_s2);
            u_i = reshape( u_n(:,end) , gridsize, gridsize);
    
            %%% (3) error analysis %%%
            measure1( j , 1) = dt;
            measure2( j , 1) = dt;
            measure1( j , 2) = norm( u_i - u_Nsq , 2) / norm(u_Nsq , 2);
            measure2( j , 2) = norm( u_i - u_Nsq , inf) / norm(u_Nsq , inf);
    
            %%% (4) computational time %%%
            measure3( j , 1) = dt;
            measure3( j , 2) = time_n;
    
            j = j + 1;
    
        end
    case 'spatial'
        
        L_s2 = utility1;
        dt = timestep;
        %%% (1) make fine grid %%%
        N_fine = gridsize(end); % number of grid points for finest grid
    
        % unit physical space domain
        x1_pts = linspace( 0 , 1 - 1/N_fine , N_fine ); 
        x2_pts = linspace( 0 , 1 - 1/N_fine , N_fine ); 
        [ x1_fine , x2_fine ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
        %%% (2) load fine solution %%%
    
        [u_n, ~, time_n] = load_2DKSsolution('time_evolution', control, dt, T, N_fine, L_s1, L_s2);
        u_Nsq = reshape( u_n(:,end) , N_fine, N_fine);
        measure1 = NaN( length(gridsize) , 2 );
        measure2 = measure1;
        measure3 = NaN( length(gridsize) , 2 );
    
        %%% (3) error analysis %%%
        measure1( end , 1) = N_fine;
        measure2( end , 1) = N_fine;
        measure1( end , 2) = 0;
        measure2( end , 2) = 0;
    
        %%% (4) computational time %%%
        measure3( end , 1) = N_fine;
        measure3( end , 2) = time_n;
    
        j = 1; % iteration counter for saved measures
    
        for i = 1 : length(gridsize)-1 % plot linf, l2 errors againt finest grid
    
            %%% (1) make coarser grid %%%
    
            N_coarse = gridsize(i); % number of grid points 
    
            % unit physical space domain
            x1_pts = linspace( 0 , 1 - 1/N_coarse , N_coarse ); 
            x2_pts = linspace( 0 , 1 - 1/N_coarse , N_coarse ); 
            [ x1_coarse , x2_coarse ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid
    
            %%% (2) load solution %%%
            [u_n, ~, time_n] = load_2DKSsolution('time_evolution', control, dt, T, N_coarse, L_s1, L_s2);
            
            u_nsq = reshape( u_n(:,end) , N_coarse, N_coarse);
            u_i = interp2( x1_coarse , x2_coarse, u_nsq, x1_fine , x2_fine); % interpolate coarse grid to compatible size
            
            %%% (3) error analysis %%%
            interpskip = ceil( length(gridsize) / i ) - 1; % interpolated grid skips this many rows/columns
            measure1( j , 1) = N_coarse;
            measure2( j , 1) = N_coarse;
            measure1( j , 2) = ...
                norm( u_i( 1:end-interpskip, 1:end-interpskip ) - u_Nsq( 1:end-interpskip, 1:end-interpskip ) , 2) ...
                / norm(u_Nsq( 1:end-interpskip, 1:end-interpskip ) , 2);
            measure2( j , 2) = ...
                norm( u_i( 1:end-interpskip, 1:end-interpskip ) - u_Nsq( 1:end-interpskip, 1:end-interpskip ) , inf) ...
                / norm(u_Nsq( 1:end-interpskip, 1:end-interpskip ) , inf);
    
            %%% (4) computational time %%%
            measure3( j , 1) = N_coarse;
            measure3( j , 2) = time_n;
    
            j = j + 1;
    
        end
    
        %{
        figure(10);
        loglog(error_2(:,1), error_2(:,2));
        title("Relative L2 Discretization Error, N = " + num2str(N_fine));
        xlabel('n = Grid Size');
        ylabel('| u_n - u_N |_2 / | u_N |_2');
    
        figure(11);
        loglog(error_inf(:,1), error_inf(:,2));
        title("Relative L-Inf Discretization Error, N = " + num2str(N_fine));
        xlabel('n = Grid Size');
        ylabel('| u_n - u_N |_inf / | u_N |_inf');
    
        figure(12);
        loglog(comptime(:,1), comptime(:,2));
        title("Computational Time");
        xlabel('Grid Size');
        ylabel('Seconds');
        %}

end
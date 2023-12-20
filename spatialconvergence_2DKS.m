function [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, method, dt, T, L_s1, L_s2)

    %%% (1) make fine grid %%%

    N_fine = gridsize(end); % number of grid points for finest grid

    % unit physical space domain
    x1_pts = linspace( 0 , 1 - 1/N_fine , N_fine ); 
    x2_pts = linspace( 0 , 1 - 1/N_fine , N_fine ); 
    [ x1_fine , x2_fine ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

    %%% (2) load fine solution %%%

    [u_n, ~, time_n] = load_2DKSsolution('time_evolution', method, dt, T, N_fine, L_s1, L_s2);
    u_Nsq = reshape( u_n(:,end) , N_fine, N_fine);
    error_2 = NaN( length(gridsize) , 2 );
    error_inf = error_2;
    comptime = NaN( length(gridsize) , 2 );

    %%% (3) error analysis %%%
    error_2( end , 1) = N_fine;
    error_inf( end , 1) = N_fine;
    error_2( end , 2) = 0;
    error_inf( end , 2) = 0;

    %%% (4) computational time %%%
    comptime( end , 1) = N_fine;
    comptime( end , 2) = time_n;

    j = 1; % iteration counter for saved measures

    for i = 1 : length(gridsize)-1 % plot linf, l2 errors againt finest grid

        %%% (1) make coarser grid %%%

        N_coarse = gridsize(i); % number of grid points 

        % unit physical space domain
        x1_pts = linspace( 0 , 1 - 1/N_coarse , N_coarse ); 
        x2_pts = linspace( 0 , 1 - 1/N_coarse , N_coarse ); 
        [ x1_coarse , x2_coarse ] = meshgrid(x1_pts,x2_pts); % 2-dimensional grid

        %%% (2) load solution %%%
        [u_n, ~, time_n] = load_2DKSsolution('time_evolution', method, dt, T, N_coarse, L_s1, L_s2);
        
        u_nsq = reshape( u_n(:,end) , N_coarse, N_coarse);
        u_i = interp2( x1_coarse , x2_coarse, u_nsq, x1_fine , x2_fine); % interpolate coarse grid to compatible size
        
        %%% (3) error analysis %%%
        interpskip = ceil( length(gridsize) / i ) - 1; % interpolated grid skips this many rows/columns
        error_2( j , 1) = N_coarse;
        error_inf( j , 1) = N_coarse;
        error_2( j , 2) = ...
            norm( u_i( 1:end-interpskip, 1:end-interpskip ) - u_Nsq( 1:end-interpskip, 1:end-interpskip ) , 2) ...
            / norm(u_Nsq( 1:end-interpskip, 1:end-interpskip ) , 2);
        error_inf( j , 2) = ...
            norm( u_i( 1:end-interpskip, 1:end-interpskip ) - u_Nsq( 1:end-interpskip, 1:end-interpskip ) , inf) ...
            / norm(u_Nsq( 1:end-interpskip, 1:end-interpskip ) , inf);

        %%% (4) computational time %%%
        comptime( j , 1) = N_coarse;
        comptime( j , 2) = time_n;

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
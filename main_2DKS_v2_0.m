%function main_2DKS_v2_0
tic

%%% choose switches %%%
run = 'IC'; % switch to 0, 'L', 'N', 'dt', 'T', 'IC', 'kappa'
test_parameter = 0; % do not alter
diagnostics = 1; % evaluating dynamics of Fourier modes
method  = 'imexrk4vec2a'; % time-stepping scheme
%IC = 'mn1'; % initial condition

%%% choose parameter ranges %%%
timewindow = linspace(0,60,7); % for time-stepping analysis
timewindow(1) = 1;
%L_scale = linspace(11,19,9)/10; % for dynamical behavior analysis
L_scale = [ 1.1 , 1.4 , 1.5 , 1.9 , 3.2 , 5.2 , 10.2]; % for dynamical behavior analysis
timestep = 10.^(-linspace(1,4,7)); % for temporal convergence analysis
gridsize = 10*4*linspace(3,21,7); % for spatial convergence analysis
initialcondition = { 'sinL' };

%%% choose default parameters %%%
%{
L_s1 = L_scale(6); % length-scale parameter in dim 1
L_s2 = L_scale(4); % length-scale parameter in dim 2
dt = timestep(5); % length of time-step
T = timewindow(2); % time window
N = gridsize(3); % number of grid points 
save_each = 100; % number of iterations between saved timepoints
%}

for lscale = 4:4

    %%% choose temporary parameters %%%
    %
    L_s1 = L_scale(lscale); % length-scale parameter in dim 1
    L_s2 = L_s1; % length-scale parameter in dim 2
    dt = 1e-2; % length of time-step
    T = 100; % time window
    N = 32; % number of grid points 
    save_each = 1/dt; % number of iterations between saved timepoints - use T*10 to save 100, 1/dt to save 1 T
    %}

    switch run % set which test to run
        case 'L'
            test_parameter = L_scale; % (dynamical behavior)
        case 'N'
            test_parameter = gridsize; % (spatial convergence)
        case 'dt'
            test_parameter = timestep; % (temporal convergence)
        case 'T'
            test_parameter = timewindow; % (temporal convergence)
        case 'IC'
            test_parameter = initialcondition; % 
        case 'kappa'
            test_parameter = 1; % (dummy variable)
            save_each = 1; % save all timesteps for adjoint solver
    end

    %
    for k = 1 : length(test_parameter) % length('x') indicates 'x' testing

        method  = 'imexrk4vec2a'; % time-stepping scheme
    
        %%% set variable parameters %%%
        switch run 
            case 'L'
                L_s2 = L_scale(k); % length-scale parameter in dim 2 
            case 'N'
                N = gridsize(k); % number of grid points 
            case 'dt'
                dt = timestep(k); % length of time-step
            case 'T'
                T = timewindow(k); % length of simulation time window
            case 'IC'
                IC = strjoin(initialcondition(k),'');
        end
        
        %%% solve PDE problem in time %%%
        %
        tic
        [ v_n , u_n ] = DirectSolve_2DKS_v2_0(IC,method,N,L_s1,L_s2,dt,T,save_each);
        time = toc;
        toc
        %}
        
        %%% compute L2 norm and Fourier mode evolution %%%
        nsave = size(u_n,2);
        normL2 = NaN(nsave,1);
        v_mean = zeros(round(sqrt((N/2)^2+(N/2)^2)) + 1,nsave);
        v_meancount = v_mean;
        for i = 1:nsave
            u_i = reshape( u_n(:,i) , [ N , N ] );
            normL2(i,1) = norm(u_i)/N;
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

        if diagnostics == 1
            % Wavenumber evolution
            wavenumberevol_file = [pwd '/data/wavenumberevol_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1) '_lY' num2str(L_s2) '.png'];
            Ntime = size(u_n,2);
            timewindow = linspace(0,T,Ntime);
            h = figure;
            semilogy(timewindow,v_mean(1,:),'LineWidth',2)
            hold on;
            for i = 2:size(v_mean,1)
                semilogy(timewindow,v_mean(i,:),'LineWidth',2)
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
            imwrite(im,wavenumberevol_file,'png');
        end

        %%% save/inspect solution %%%
        switch run 
            case {'L','N','dt','T','IC'}                
                %
                method = IC;
                save_2DKSsolution('time_evolution', u_n, v_n, time, method, dt, T, N, L_s1, L_s2); % save solution
                %{
                [u_n, v_n, ~] = load_2DKSsolution('time_evolution', IC, dt, T, N, L_s1, L_s2); % load solution
                %}
                %plot2DKS(v_n , u_n, 'initial', method, N, dt, T, L_s1, L_s2); % save/inspect initial state
                plot2DKS(v_n , u_n, 'terminal', method, N, dt, T, L_s1, L_s2, normL2, v_mean); % save/inspect terminal state
                plot2DKS(v_n , u_n, 'gif', method, N, dt, T, L_s1, L_s2, normL2, v_mean); % save/inspect surface time evolution
                close all
                %
        end

    end
    %

    %{
    switch run 
        case 'N' % spatial convergence: error analysis and computational time
            [error_2,error_inf,comptime] = spatialconvergence_2DKS(gridsize, method, dt, T, L_s1, L_s2);
            save_measures('spatial', error_2, error_inf, comptime, method, 0, dt, T, L_s1, L_s2);
        case 'dt' % temporal convergence: error analysis and computational time
            [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, method, N, T, L_s1, L_s2);
            save_measures('temporal', error_2, error_inf, comptime, method, N, 0, T, L_s1, L_s2);
    end
    %}

end

%{
dt = 1e-3; % length of time-step
T = 150; % time window
N = 512; % number of grid points 
save_each = 1/dt;

method  = 'imexrk4vec2a'; % time-stepping scheme
IC = 'mn5'; % initial condition
%L_s1 = 3.2; % length-scale parameter in dim 1
L_s2 = L_s1; % length-scale parameter in dim 2
tic
[ v_n , u_n ] = DirectSolve_2DKS_v2_0(IC,method,N,L_s1,L_s2,dt,T,save_each);
time = toc;
toc
method = IC;
save_2DKSsolution('time_evolution', u_n, v_n, time, method, dt, T, N, L_s1, L_s2); % save solution
plot2DKS(v_n , u_n, 'terminal', method, N, dt, T, L_s1, L_s2); % save/inspect terminal state
plot2DKS(v_n , u_n, 'gif', method, N, dt, T, L_s1, L_s2); % save/inspect surface time evolution
close all

%}
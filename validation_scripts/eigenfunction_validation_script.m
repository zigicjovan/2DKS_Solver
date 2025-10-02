close all
L_scale = sqrt(32);                                % domain sizes
%[sqrt(2),2,sqrt(8),sqrt(10),4,sqrt(18),sqrt(20),sqrt(26),sqrt(32)]
%[sqrt(3),sqrt(6),3,sqrt(13),sqrt(17),sqrt(19),sqrt(23),sqrt(29)]
timestep = .005;                                % time-step sizes
gridsize = 48;                                  % grid sizes
timewindow = 30;                                % time windows
initialcondition = {'noise4'};                      % initial conditions
tol = 1e-10;                                     % set optimization tolerance critera
K = 0;

%%% choose default parameters %%%
L_s1 = L_scale(1);                              % length-scale parameter in dim 1
L_s2 = L_scale(1);                              % length-scale parameter in dim 2
dt = timestep(1);                               % length of time-step
N = gridsize(1);                                % number of grid points
T = timewindow(1);                              % length of simulation time window
IC = strjoin(initialcondition(1),'');           % initial condition

u_IC_opt = load_2DKSsolution('optimal', IC, dt, T, N, K, L_s1, L_s2, tol, 0);
[match_scored,ampstarsd,modesd] = validation_script(u_IC_opt,L_s1, N, T,IC,'dominant');
%[match_scorea,ampstarsa,modesa] = validation_script(u_IC_opt,L_s1, N, T,IC,'active');
%[match_scoref,ampstarsf,modesf] = validation_script(u_IC_opt,L_s1, N, T,IC,'full');
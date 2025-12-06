%% RCG v RF eigenfunction validation

rg141 = load([pwd '/data/optimization/diagnostics_noise_optimized_N_16_dt_0.001_K_0_Ls1_1.41_Ls2_1.41_T_10_tol_1e-10_rg.dat']);
rcg141 = load([pwd '/data/optimization/diagnostics_noise_optimized_N_16_dt_0.001_K_0_Ls1_1.41_Ls2_1.41_T_10_tol_1e-10_rcg.dat']);

h = figure;
set(gcf,'color','white')
set(gca,'color','white')  
sgtitle('RCG vs. RG algorithm convergence for $K=10^{-4},\ell=\sqrt{2},T=10$','Interpreter','latex')
K = 1e-4;
T = 10;
ell = sqrt(2);
k1 = 1; k2 = 0;
lam = (k1/ell)^2*(1-(k1/ell)^2) + (k2/ell)^2*(1-(k2/ell)^2) - 2*((k1*k2)/(ell*ell))^2;

subplot(3,2,1:2)
semilogy(rg141(:,1),'r-x')
hold on
semilogy(rcg141(:,1),'b-o')
yline((K*exp(2*lam*T)),'--','LineWidth',1)
ylim([rg141(1,1) 1.5*rcg141(end,1)])
hold off
legendlist(1) = {'RG'};
legendlist(2) = {'RCG'};
legendlist(3) = {'$Ke^{2\lambda^* T}$'};
legend(legendlist,'Interpreter','latex','Location','southeast','Orientation','horizontal','Box','off','FontSize',14)
xlabel('Iteration $n$','Interpreter','latex','FontSize',14); 
ylabel('$\| \phi(T;{\varphi^{(n)}}) \|^2_{L^2}$','Interpreter','latex','FontSize',14);

subplot(3,2,3:4)
semilogx(rg141(:,1),'r-x')
hold on
semilogx(rcg141(:,1),'b-o')
yline((K*exp(2*lam*T)),'--','LineWidth',1)
ylim([rcg141(3,1)-1e-5 (K*exp(2*lam*T)+1e-5)])
hold off
legend(legendlist,'Interpreter','latex','Location','southeast','Orientation','horizontal','Box','off','FontSize',14)
xlabel('Iteration $n$','Interpreter','latex','FontSize',14); 
ylabel('$\| \phi(T;{\varphi^{(n)}}) \|^2_{L^2}$','Interpreter','latex','FontSize',14);

subplot(3,2,5:6)
semilogx(rg141(:,1),'r-x')
hold on
semilogx(rcg141(:,1),'b-o')
yline((K*exp(2*lam*T)),'--','LineWidth',1)
ylim([rcg141(3,1) rcg141(end,1)+1e-9])
hold off
legend(legendlist,'Interpreter','latex','Location','southeast','Orientation','horizontal','Box','off','FontSize',14)
xlabel('Iteration $n$','Interpreter','latex','FontSize',14); 
ylabel('$\| \phi(T;{\varphi^{(n)}}) \|^2_{L^2}$','Interpreter','latex','FontSize',14);

%{
subplot(2,2,4)
semilogx(rg141(:,1),'r-x')
hold on
semilogx(rcg141(:,1),'b-o')
%ylim([rcg141(6,1)-1e-9 rg141(end,1)+1e-11])
%ylim([rcg141(6,1)-1e-9 (K*exp(2*lam*T))+1e-11])
hold off
legendlist(1) = {'RG'};
legendlist(2) = {'RCG'};
legend(legendlist,'Interpreter','latex','Location','northeast','NumColumns',1,'Box','off','FontSize',14)
xlabel('Iteration $n$','Interpreter','latex','FontSize',14); 
%ylabel('$\| \phi(T;{\varphi^{(n)}}) \|^2_{L^2}$','Interpreter','latex','FontSize',14);
%}
filename = [pwd '/media/optimization/RCGvRG_lmax'];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

rcg141 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_1.41_Ls2_1.41_T_30_tol_1e-10.dat']);
rcg200 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_2.00_Ls2_2.00_T_30_tol_1e-10.dat']);
rcg283 =  load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_2.83_Ls2_2.83_T_30_tol_1e-10.dat']);
rcg316 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_3.16_Ls2_3.16_T_30_tol_1e-10.dat']);
rcg400 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_4.00_Ls2_4.00_T_30_tol_1e-10.dat']);
rcg424 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_4.24_Ls2_4.24_T_30_tol_1e-10.dat']);
%rcg436 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_4.36_Ls2_4.36_T_30_tol_1e-10.dat']);
rcg447 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_4.47_Ls2_4.47_T_30_tol_1e-10.dat']);
rcg510 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_5.10_Ls2_5.10_T_30_tol_1e-10.dat']);
rcg566 = load([pwd '/data/optimization/ef_validation/diagnostics_noise4_optimized_N_48_dt_0.005_K_0_Ls1_5.66_Ls2_5.66_T_30_tol_1e-10.dat']);
rg141 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_1.414_Ls2_1.414_tol_1e-10.dat']);
rg200 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_2.000_Ls2_2.000_tol_1e-10.dat']);
rg283 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_2.828_Ls2_2.828_tol_1e-10.dat']);
rg316 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_3.162_Ls2_3.162_tol_1e-10.dat']);
rg400 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_4.000_Ls2_4.000_tol_1e-10.dat']);
rg424 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_4.243_Ls2_4.243_tol_1e-10.dat']);
%rg436 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_4.359_Ls2_4.359_tol_1e-10.dat']);
rg447 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_4.472_Ls2_4.472_tol_1e-10.dat']);
rg510 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_5.099_Ls2_5.099_tol_1e-10.dat']);
rg566 = load([pwd '/data/optimization/ef_validation/diagnostics_optimized_noise4_N_48_T_30_dt_0.005_Ls1_5.657_Ls2_5.657_tol_1e-10.dat']);

rcgJ = NaN(1001,9);
rgJ = rcgJ;

rcgJ(1:size(rcg141,1),1) = rcg141(:,1);
rcgJ(1:size(rcg200,1),2) = rcg200(:,1);
rcgJ(1:size(rcg283,1),3) = rcg283(:,1);
rcgJ(1:size(rcg316,1),4) = rcg316(:,1);
rcgJ(1:size(rcg400,1),5) = rcg400(:,1);
rcgJ(1:size(rcg424,1),6) = rcg424(:,1);
rcgJ(1:size(rcg447,1),7) = rcg447(:,1);
rcgJ(1:size(rcg510,1),8) = rcg510(:,1);
rcgJ(1:size(rcg566,1),9) = rcg566(:,1);

rgJ(1:size(rg141,1),1) = rg141(:,1);
rgJ(1:size(rg200,1),2) = rg200(:,1);
rgJ(1:size(rg283,1),3) = rg283(:,1);
rgJ(1:size(rg316,1),4) = rg316(:,1);
rgJ(1:size(rg400,1),5) = rg400(:,1);
rgJ(1:size(rg424,1),6) = rg424(:,1);
rgJ(1:size(rg447,1),7) = rg447(:,1);
rgJ(1:size(rg510,1),8) = rg510(:,1);
rgJ(1:size(rg566,1),9) = rg566(:,1);


h = figure;
plot(rgJ,'LineWidth',0.5,'Marker','.','Color','b')
hold on
plot(rcgJ,'LineWidth',0.5,'Marker','.','Color','r')
hold off
set(gcf,'Position',[100 100 900 750])
xlabel('Iteration number','Interpreter','latex'); 
ylabel('$\| {\phi(T;\tilde\varphi)} \|_{L^2}$','Interpreter','latex');
fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title('RCG (red) vs. RG (blue) algorithm convergence','Interpreter','latex')
%subtitle(optparfiglist,'Interpreter','latex','FontSize',14)
%legend(['$\phi(t,\varphi_{' originalIC '})$'],'$\phi(t,\widetilde{\varphi})$','Interpreter','latex','Location','northwest')
filename = [pwd '/media/optimization/RCGvRG_lmax2'];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% Extend branches by combining saved data files with current data file

Jinitdata1 = load([pwd '/data/optimization/Jinit_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);
Joptdata1 = load([pwd '/data/optimization/Jopt_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);

Jinitdatanew2 = [ Jinitdatanew(1:22,:) ; Jinitdata(1:2,:) ; ...
    Jinitdatanew(23:44,:) ; Jinitdata(3:4,:) ; ...
    Jinitdatanew(45:66,:) ; Jinitdata(5:6,:) ; ...
    Jinitdatanew(67:88,:) ; Jinitdata(7:8,:) ; ...
    Jinitdatanew(89:110,:) ; Jinitdata(9:10,:) ];

Joptdatanew2 = [ Joptdatanew(1:22,:) ; Joptdata(1:2,:) ; ...
    Joptdatanew(23:44,:) ; Joptdata(3:4,:) ; ...
    Joptdatanew(45:66,:) ; Joptdata(5:6,:) ; ...
    Joptdatanew(67:88,:) ; Joptdata(7:8,:) ; ...
    Joptdatanew(89:110,:) ; Joptdata(9:10,:) ];

plot_measures('optimization', dt, initialcondition, N, T, K, L_s1, L_s2, Jinitdatanew2, Joptdatanew2, IC, tol);
save_measures('optimization', Jinitdatanew2, Joptdatanew2, numberoftests, initialcondition, N, dt, timewindow, K, L_scale, L_s2);

%% compare to s1 benchmark
%Jinitdata_s1 = load([pwd '/data/optimization/Jinit_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);
Joptdata_s1 = load([pwd '/data/optimization/Jopt_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);

%Jinitdata_stg1 = load([pwd '/data/optimization/Jinit_stg1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);
Joptdata_stg1 = load([pwd '/data/optimization/Jopt_stg1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);

Joptdata_noise = load([pwd '/data/optimization/Jopt_noise_N_32_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_150.00.dat']);

Joptdata_stg30 = load([pwd '/data/optimization/Jopt_stg30_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_5.07_150.00.dat']);

h = figure;

plot(Joptdata_s1(14:24,3),Joptdata_s1(14:24,4),'b-','MarkerSize',5)
hold on
plot(Joptdata_stg1(14:24,3),Joptdata_stg1(14:24,4),'r--','MarkerSize',5)
plot(Joptdata_stg30(2:12,3),Joptdata_stg30(2:12,4),'kx','MarkerSize',7)
plot(Joptdata_noise(14:24,3),Joptdata_noise(14:24,4),'go','MarkerSize',4)

plot(Joptdata_s1(38:48,3),Joptdata_s1(38:48,4),'b-','MarkerSize',5)
plot(Joptdata_s1(62:72,3),Joptdata_s1(62:72,4),'b-','MarkerSize',5)
plot(Joptdata_s1(86:96,3),Joptdata_s1(86:96,4),'b-','MarkerSize',5)
plot(Joptdata_s1(110:120,3),Joptdata_s1(110:120,4),'b-','MarkerSize',5)

plot(Joptdata_stg1(38:48,3),Joptdata_stg1(38:48,4),'r--','MarkerSize',5)
plot(Joptdata_stg1(62:72,3),Joptdata_stg1(62:72,4),'r--','MarkerSize',5)
plot(Joptdata_stg1(86:96,3),Joptdata_stg1(86:96,4),'r--','MarkerSize',5)
plot(Joptdata_stg1(110:120,3),Joptdata_stg1(110:120,4),'r--','MarkerSize',5)

plot(Joptdata_stg30(14:24,3),Joptdata_stg30(14:24,4),'kx','MarkerSize',7)
plot(Joptdata_stg30(26:36,3),Joptdata_stg30(26:36,4),'kx','MarkerSize',7)
plot(Joptdata_stg30(38:48,3),Joptdata_stg30(38:48,4),'kx','MarkerSize',7)
plot(Joptdata_stg30(50:60,3),Joptdata_stg30(50:60,4),'kx','MarkerSize',7)

plot(Joptdata_noise(38:48,3),Joptdata_noise(38:48,4),'go','MarkerSize',4)
plot(Joptdata_noise(62:72,3),Joptdata_noise(62:72,4),'go','MarkerSize',4)
plot(Joptdata_noise(86:96,3),Joptdata_noise(86:96,4),'go','MarkerSize',4)
plot(Joptdata_noise(110:120,3),Joptdata_noise(110:120,4),'go','MarkerSize',4)

hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1200 900])
xlabel('Time Window $T$','Interpreter','latex'); 
ylabel('$\| {\phi(T;\varphi)} \|^2_{L^2}$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = '$N = 48, {\Delta}t = 0.005, K = 1, \ell = [1.02, 1.10]$';
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
legendlist(1) = {['$\| \phi(T;\tilde{\varphi}_{' num2str(1,'%.0f') ',2\pi(\ell),T}) \|^2_{L^2}~,~\varphi^{(0)}=\varphi_{s1}~~~~~~~~$']};
legendlist(2) = {['$\| \phi(T;\tilde{\varphi}_{' num2str(1,'%.0f') ',2\pi(\ell),T}) \|^2_{L^2}~,~\varphi^{(0)}=\varphi_{stg1}~~~~~~~~$']};
legendlist(3) = {['$\| \phi(T;\tilde{\varphi}_{' num2str(1,'%.0f') ',2\pi(\ell),T}) \|^2_{L^2}~,~\varphi^{(0)}=\varphi_{stg30}~~~~~~~~$']};
legendlist(4) = {['$\| \phi(T;\tilde{\varphi}_{' num2str(1,'%.0f') ',2\pi(\ell),T}) \|^2_{L^2}~,~\varphi^{(0)}=\varphi_{noise}~~~~~~~~$']};
legend(legendlist,'Interpreter','latex','Location','southoutside','NumColumns',2,'Box','off')
hold off
filename = [pwd '/media/optimization/branches_comp_4ic' ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])


%% max T K0-K5, 1.02-1.10
Joptdata_K0 = load([pwd '/data/optimization/Jopt_s1_N_24_dt_0.005_K_1_1_L_1.02_1.10_T_1.00000_100.00000.dat']);
Joptdata_K1 = load([pwd '/data/optimization/Jopt_s1_N_32_dt_0.0025_K_10_10_L_1.02_1.10_T_0.31623_31.62278.dat']);
Joptdata_K2 = load([pwd '/data/optimization/Jopt_s1_N_32_dt_0.0025_K_100_100_L_1.02_1.10_T_0.10000_10.00000.dat']);
Joptdata_K3 = load([pwd '/data/optimization/Jopt_s1_N_48_dt_0.001_K_1000_1000_L_1.02_1.10_T_0.03162_3.16228.dat']);
Joptdata_K4 = load([pwd '/data/optimization/Jopt_s1_N_64_dt_0.0001_K_10000_10000_L_1.02_1.10_T_0.01000_0.10000.dat']);

h = figure;

loglog(Joptdata_K0(1:21,3),Joptdata_K0(1:21,4),'r-o','MarkerSize',5)
hold on
loglog(Joptdata_K0(22:42,3),Joptdata_K0(22:42,4),'b-o','MarkerSize',5)
loglog(Joptdata_K0(43:63,3),Joptdata_K0(43:63,4),'g-o','MarkerSize',5)
loglog(Joptdata_K0(64:84,3),Joptdata_K0(64:84,4),'k-o','MarkerSize',5)
loglog(Joptdata_K0(85:105,3),Joptdata_K0(85:105,4),'m-o','MarkerSize',5)

loglog(Joptdata_K1(1:21,3),Joptdata_K1(1:21,4),'r-o','MarkerSize',5)
loglog(Joptdata_K1(22:42,3),Joptdata_K1(22:42,4),'b-o','MarkerSize',5)
loglog(Joptdata_K1(43:63,3),Joptdata_K1(43:63,4),'g-o','MarkerSize',5)
loglog(Joptdata_K1(64:84,3),Joptdata_K1(64:84,4),'k-o','MarkerSize',5)
loglog(Joptdata_K1(85:105,3),Joptdata_K1(85:105,4),'m-o','MarkerSize',5)

loglog(Joptdata_K2(1:21,3),Joptdata_K2(1:21,4),'r-o','MarkerSize',5)
loglog(Joptdata_K2(22:42,3),Joptdata_K2(22:42,4),'b-o','MarkerSize',5)
loglog(Joptdata_K2(43:63,3),Joptdata_K2(43:63,4),'g-o','MarkerSize',5)
loglog(Joptdata_K2(64:84,3),Joptdata_K2(64:84,4),'k-o','MarkerSize',5)
loglog(Joptdata_K2(85:105,3),Joptdata_K2(85:105,4),'m-o','MarkerSize',5)

loglog(Joptdata_K3(1:21,3),Joptdata_K3(1:21,4),'r-o','MarkerSize',5)
loglog(Joptdata_K3(22:42,3),Joptdata_K3(22:42,4),'b-o','MarkerSize',5)
loglog(Joptdata_K3(43:63,3),Joptdata_K3(43:63,4),'g-o','MarkerSize',5)
loglog(Joptdata_K3(64:84,3),Joptdata_K3(64:84,4),'k-o','MarkerSize',5)
loglog(Joptdata_K3(85:105,3),Joptdata_K3(85:105,4),'m-o','MarkerSize',5)

loglog(Joptdata_K4(1:11,3),Joptdata_K4(1:11,4),'r-o','MarkerSize',5)
loglog(Joptdata_K4(12:22,3),Joptdata_K4(12:22,4),'b-o','MarkerSize',5)
loglog(Joptdata_K4(23:33,3),Joptdata_K4(23:33,4),'g-o','MarkerSize',5)
loglog(Joptdata_K4(34:44,3),Joptdata_K4(34:44,4),'k-o','MarkerSize',5)
loglog(Joptdata_K4(45:55,3),Joptdata_K4(45:55,4),'m-o','MarkerSize',5)

hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1200 900])
xlabel('Time Window $T$','Interpreter','latex'); 
ylabel('$\| {\phi(T;\varphi)} \|^2_{L^2}$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = '$\varphi^{(0)}=\varphi_{s1}, K = \left[10^{0}, 10^{4}\right], \ell = [1.02, 1.10]$';
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
legendlist(1) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.02),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(2) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.04),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(3) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.06),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(4) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.08),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(5) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.10),T}) \|^2_{L^2}~~~~~~~~$'};
legend(legendlist,'Interpreter','latex','Location','northeast','NumColumns',1,'Box','off')

filename = [pwd '/media/optimization/branches_comp_5K' ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% max T K4-K6, 1.02-1.40
ellarray = 1.02:.02:1.40;
numK = 5;
Joptdata_K40 = NaN(30,length(ellarray)+1);
Joptdata_K45 = NaN(30,length(ellarray)+1);
Joptdata_K50 = NaN(30,length(ellarray)+1);
Joptdata_K55 = NaN(30,length(ellarray)+1);
Joptdata_K60 = NaN(30,length(ellarray)+1);
Jmaxopt = NaN(length(ellarray)*numK,3);
Jmaxoptplot = NaN(numK,length(ellarray)+1);
for i = 1:length(ellarray)
    Joptdata_K40temp = load([pwd '/data/branches/Jopt_s1_N_96_dt_0.0001_K_10000_10000_L_' num2str(ellarray(i),'%.2f') '_' num2str(ellarray(i),'%.2f') '.dat']);
    Joptdata_K45temp = load([pwd '/data/branches/Jopt_s1_N_144_dt_2e-05_K_31623_31623_L_' num2str(ellarray(i),'%.2f') '_' num2str(ellarray(i),'%.2f') '.dat']);
    Joptdata_K50temp = load([pwd '/data/branches/Jopt_s1_N_256_dt_2e-06_K_100000_100000_L_' num2str(ellarray(i),'%.2f') '_' num2str(ellarray(i),'%.2f') '.dat']);
    Joptdata_K55temp = load([pwd '/data/branches/Jopt_s1_N_320_dt_5e-07_K_316228_316228_L_' num2str(ellarray(i),'%.2f') '_' num2str(ellarray(i),'%.2f') '.dat']);
    Joptdata_K60temp = load([pwd '/data/branches/Jopt_s1_N_360_dt_2e-07_K_1000000_1000000_L_' num2str(ellarray(i),'%.2f') '_' num2str(ellarray(i),'%.2f') '.dat']);
    Joptdata_K40(1:size(Joptdata_K40temp,1),i+1) = Joptdata_K40temp(:,4);    
    Joptdata_K45(1:size(Joptdata_K45temp,1),i+1) = Joptdata_K45temp(:,4);
    Joptdata_K50(1:size(Joptdata_K50temp,1),i+1) = Joptdata_K50temp(:,4);
    Joptdata_K55(1:size(Joptdata_K55temp,1),i+1) = Joptdata_K55temp(:,4);
    Joptdata_K60(1:size(Joptdata_K60temp,1),i+1) = Joptdata_K60temp(:,4);
    Jmaxopt((i-1)*numK+1,:) = [ Joptdata_K40temp(1,1:2), max(Joptdata_K40temp(:,4)) ];
    Jmaxopt((i-1)*numK+2,:) = [ Joptdata_K45temp(1,1:2), max(Joptdata_K45temp(:,4)) ];
    Jmaxopt((i-1)*numK+3,:) = [ Joptdata_K50temp(1,1:2), max(Joptdata_K50temp(:,4)) ];
    Jmaxopt((i-1)*numK+4,:) = [ Joptdata_K55temp(1,1:2), max(Joptdata_K55temp(:,4)) ];
    Jmaxopt((i-1)*numK+5,:) = [ Joptdata_K60temp(1,1:2), max(Joptdata_K60temp(:,4)) ];
    Jmaxoptplot(:,i+1) = Jmaxopt((i-1)*numK+1:(i-1)*numK+5,3);
end
Joptdata_K40(1:size(Joptdata_K40temp,1),1) = Joptdata_K40temp(:,3);
Joptdata_K45(1:size(Joptdata_K45temp,1),1) = Joptdata_K45temp(:,3);
Joptdata_K50(1:size(Joptdata_K50temp,1),1) = Joptdata_K50temp(:,3);
Joptdata_K55(1:size(Joptdata_K55temp,1),1) = Joptdata_K55temp(:,3);
Joptdata_K60(1:size(Joptdata_K60temp,1),1) = Joptdata_K60temp(:,3);
Jmaxoptplot(:,1) = Jmaxopt(1:numK,1);

h = figure;

loglog(Joptdata_K40(:,1),Joptdata_K40(:,2:length(ellarray)),'-o','MarkerSize',5)
hold on
loglog(Joptdata_K45(:,1),Joptdata_K45(:,2:length(ellarray)),'-o','MarkerSize',5)
loglog(Joptdata_K50(:,1),Joptdata_K50(:,2:length(ellarray)),'-o','MarkerSize',5)
loglog(Joptdata_K55(:,1),Joptdata_K55(:,2:length(ellarray)),'-o','MarkerSize',5)
loglog(Joptdata_K60(:,1),Joptdata_K60(:,2:length(ellarray)),'-o','MarkerSize',5)

hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1250 700])
xlabel('Time Window $T$','Interpreter','latex'); 
ylabel('$\| {\phi(T;\varphi)} \|^2_{L^2}$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = '$K = \left[10^{4},10^{4.5}, ... , 10^{6}\right], \ell = [1.02,1.04, ... , 1.40]$';
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
%{
legendlist(1) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.02),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(2) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.04),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(3) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.06),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(4) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.08),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(5) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.10),T}) \|^2_{L^2}~~~~~~~~$'};
legend(legendlist,'Interpreter','latex','Location','northeast','NumColumns',1,'Box','off')
%}
filename = [pwd '/media/optimization/branches_K4-6' ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

h = figure;

xlin = logspace(4,6,5);
loglog(xlin, 10^(-2.05).*(xlin.^(1.9/1.0)),'--','MarkerSize',5)
hold on
loglog(Jmaxoptplot(:,1),Jmaxoptplot(:,2:length(ellarray)),'-o','MarkerSize',5)
hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1250 700])
xlabel('Initial Energy $K$','Interpreter','latex'); 
ylabel('Peak Terminal Energy $\tilde K$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = '$K = \left[10^{4},10^{4.5}, ... , 10^{6}\right], \ell = [1.02,1.04, ... , 1.40]$';
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
legend('$K^{\alpha},\alpha = 1.9$~','Interpreter','latex','Location','southeast','Box','off','FontSize',18)
%{
legendlist(1) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.02),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(2) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.04),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(3) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.06),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(4) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.08),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(5) = {'$\| \phi(T;\tilde{\varphi}_{K,2\pi(1.10),T}) \|^2_{L^2}~~~~~~~~$'};
legend(legendlist,'Interpreter','latex','Location','northeast','NumColumns',1,'Box','off')
%}
filename = [pwd '/media/optimization/powerlaw_K4-6' ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% long T K0-K5, 1.02-1.10

en102_0 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_1_Ls1_1.02_Ls2_1.02_T_100_tol_1e-08_RCG.dat']);
en102_1 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_10_Ls1_1.02_Ls2_1.02_T_100_tol_1e-10_RCG.dat']);
en102_2 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_100_Ls1_1.02_Ls2_1.02_T_100_tol_1e-10_RCG.dat']);
en102_3 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_48_dt_0.0005_K_1000_Ls1_1.02_Ls2_1.02_T_50_tol_1e-10_RCG.dat']);
en102_4 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_64_dt_5e-05_K_10000_Ls1_1.02_Ls2_1.02_T_5_tol_1e-08_RCG.dat']);

en104_0 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_1_Ls1_1.04_Ls2_1.04_T_100_tol_1e-08_RCG.dat']);
en104_1 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_10_Ls1_1.04_Ls2_1.04_T_100_tol_1e-10_RCG.dat']);
en104_2 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_100_Ls1_1.04_Ls2_1.04_T_100_tol_1e-10_RCG.dat']);
en104_3 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_48_dt_0.0005_K_1000_Ls1_1.04_Ls2_1.04_T_50_tol_1e-10_RCG.dat']);
en104_4 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_64_dt_5e-05_K_10000_Ls1_1.04_Ls2_1.04_T_5_tol_1e-08_RCG.dat']);

en106_0 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.002_K_1_Ls1_1.06_Ls2_1.06_T_100_tol_1e-08_RCG.dat']);
en106_1 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_10_Ls1_1.06_Ls2_1.06_T_100_tol_1e-10_RCG.dat']);
en106_2 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_100_Ls1_1.06_Ls2_1.06_T_100_tol_1e-10_RCG.dat']);
en106_3 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_48_dt_0.0005_K_1000_Ls1_1.06_Ls2_1.06_T_50_tol_1e-10_RCG.dat']);
en106_4 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_64_dt_5e-05_K_10000_Ls1_1.06_Ls2_1.06_T_5_tol_1e-08_RCG.dat']);

en108_0 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.002_K_1_Ls1_1.08_Ls2_1.08_T_100_tol_1e-08_RCG.dat']);
en108_1 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_10_Ls1_1.08_Ls2_1.08_T_100_tol_1e-10_RCG.dat']);
en108_2 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_100_Ls1_1.08_Ls2_1.08_T_100_tol_1e-10_RCG.dat']);
en108_3 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_48_dt_0.0005_K_1000_Ls1_1.08_Ls2_1.08_T_50_tol_1e-10_RCG.dat']);
en108_4 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_64_dt_5e-05_K_10000_Ls1_1.08_Ls2_1.08_T_5_tol_1e-08_RCG.dat']);

en110_0 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0015_K_1_Ls1_1.10_Ls2_1.10_T_100_tol_1e-08_RCG.dat']);
en110_1 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_10_Ls1_1.10_Ls2_1.10_T_100_tol_1e-10_RCG.dat']);
en110_2 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_32_dt_0.0025_K_100_Ls1_1.10_Ls2_1.10_T_100_tol_1e-10_RCG.dat']);
en110_3 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_48_dt_0.0005_K_1000_Ls1_1.10_Ls2_1.10_T_50_tol_1e-10_RCG.dat']);
en110_4 = load([pwd '/data/energyL2/energyL2_s1_optimized_N_64_dt_5e-05_K_10000_Ls1_1.10_Ls2_1.10_T_5_tol_1e-08_RCG.dat']);

h = figure;
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[0 0 1250 700])
sgtitle('Asymptotic Behaviour of $L^2$ Energy-Optimized 2D Kuramoto-Sivashinsky Solutions','Interpreter','latex','FontSize',18)

subplot(1,2,1)
semilogy(linspace(0,5,size(en102_4,1)),en102_4,'r-','MarkerSize',5)
hold on
semilogy(linspace(0,5,size(en104_4,1)),en104_4,'b-','MarkerSize',5)
semilogy(linspace(0,5,size(en106_4,1)),en106_4,'g-','MarkerSize',5)
semilogy(linspace(0,5,size(en108_4,1)),en108_4,'k-','MarkerSize',5)
semilogy(linspace(0,5,size(en110_4,1)),en110_4,'m-','MarkerSize',5)

semilogy(linspace(0,50,size(en102_3,1)),en102_3,'r-','MarkerSize',5)
hold on
semilogy(linspace(0,50,size(en104_3,1)),en104_3,'b-','MarkerSize',5)
semilogy(linspace(0,50,size(en106_3,1)),en106_3,'g-','MarkerSize',5)
semilogy(linspace(0,50,size(en108_3,1)),en108_3,'k-','MarkerSize',5)
semilogy(linspace(0,50,size(en110_3,1)),en110_3,'m-','MarkerSize',5)

semilogy(linspace(0,100,size(en102_2,1)),en102_2,'r-','MarkerSize',5)
hold on
semilogy(linspace(0,100,size(en104_2,1)),en104_2,'b-','MarkerSize',5)
semilogy(linspace(0,100,size(en106_2,1)),en106_2,'g-','MarkerSize',5)
semilogy(linspace(0,100,size(en108_2,1)),en108_2,'k-','MarkerSize',5)
semilogy(linspace(0,100,size(en110_2,1)),en110_2,'m-','MarkerSize',5)

semilogy(linspace(0,100,size(en102_1,1)),en102_1,'r-','MarkerSize',5)
hold on
semilogy(linspace(0,100,size(en104_1,1)),en104_1,'b-','MarkerSize',5)
semilogy(linspace(0,100,size(en106_1,1)),en106_1,'g-','MarkerSize',5)
semilogy(linspace(0,100,size(en108_1,1)),en108_1,'k-','MarkerSize',5)
semilogy(linspace(0,100,size(en110_1,1)),en110_1,'m-','MarkerSize',5)

semilogy(linspace(0,100,size(en102_0,1)),en102_0,'r-','MarkerSize',5)
hold on
semilogy(linspace(0,100,size(en104_0,1)),en104_0,'b-','MarkerSize',5)
semilogy(linspace(0,100,size(en106_0,1)),en106_0,'g-','MarkerSize',5)
semilogy(linspace(0,100,size(en108_0,1)),en108_0,'k-','MarkerSize',5)
semilogy(linspace(0,100,size(en110_0,1)),en110_0,'m-','MarkerSize',5)

xlabel('Time Window $T$','Interpreter','latex','FontSize',16); 
ylabel('$\| {\phi(t;\varphi)} \|^2_{L^2}$','Interpreter','latex','FontSize',16);

%parfiglistInterval = '$\varphi^{(0)}=\varphi_{s1}, K = \left[10^{0}, 10^{4}\right], \ell = [1.02, 1.10]$';
%subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
legendlist(1) = {'$\| \phi(t;\max_T \tilde{\varphi}_{K,2\pi(1.02),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(2) = {'$\| \phi(t;\max_T \tilde{\varphi}_{K,2\pi(1.04),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(3) = {'$\| \phi(t;\max_T \tilde{\varphi}_{K,2\pi(1.06),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(4) = {'$\| \phi(t;\max_T \tilde{\varphi}_{K,2\pi(1.08),T}) \|^2_{L^2}~~~~~~~~$'};
legendlist(5) = {'$\| \phi(t;\max_T \tilde{\varphi}_{K,2\pi(1.10),T}) \|^2_{L^2}~~~~~~~~$'};
legend(legendlist,'Interpreter','latex','Location','northeast','NumColumns',1,'Box','off','FontSize',16)
set(gca,'FontSize',16)

subplot(1,2,2)
loglog(linspace(0,5,size(en102_4,1)),en102_4,'r-','MarkerSize',5)
hold on
loglog(linspace(0,5,size(en104_4,1)),en104_4,'b-','MarkerSize',5)
loglog(linspace(0,5,size(en106_4,1)),en106_4,'g-','MarkerSize',5)
loglog(linspace(0,5,size(en108_4,1)),en108_4,'k-','MarkerSize',5)
loglog(linspace(0,5,size(en110_4,1)),en110_4,'m-','MarkerSize',5)

loglog(linspace(0,50,size(en102_3,1)),en102_3,'r-','MarkerSize',5)
hold on
loglog(linspace(0,50,size(en104_3,1)),en104_3,'b-','MarkerSize',5)
loglog(linspace(0,50,size(en106_3,1)),en106_3,'g-','MarkerSize',5)
loglog(linspace(0,50,size(en108_3,1)),en108_3,'k-','MarkerSize',5)
loglog(linspace(0,50,size(en110_3,1)),en110_3,'m-','MarkerSize',5)

loglog(linspace(0,100,size(en102_2,1)),en102_2,'r-','MarkerSize',5)
hold on
loglog(linspace(0,100,size(en104_2,1)),en104_2,'b-','MarkerSize',5)
loglog(linspace(0,100,size(en106_2,1)),en106_2,'g-','MarkerSize',5)
loglog(linspace(0,100,size(en108_2,1)),en108_2,'k-','MarkerSize',5)
loglog(linspace(0,100,size(en110_2,1)),en110_2,'m-','MarkerSize',5)

loglog(linspace(0,100,size(en102_1,1)),en102_1,'r-','MarkerSize',5)
hold on
loglog(linspace(0,100,size(en104_1,1)),en104_1,'b-','MarkerSize',5)
loglog(linspace(0,100,size(en106_1,1)),en106_1,'g-','MarkerSize',5)
loglog(linspace(0,100,size(en108_1,1)),en108_1,'k-','MarkerSize',5)
loglog(linspace(0,100,size(en110_1,1)),en110_1,'m-','MarkerSize',5)

loglog(linspace(0,100,size(en102_0,1)),en102_0,'r-','MarkerSize',5)
hold on
loglog(linspace(0,100,size(en104_0,1)),en104_0,'b-','MarkerSize',5)
loglog(linspace(0,100,size(en106_0,1)),en106_0,'g-','MarkerSize',5)
loglog(linspace(0,100,size(en108_0,1)),en108_0,'k-','MarkerSize',5)
loglog(linspace(0,100,size(en110_0,1)),en110_0,'m-','MarkerSize',5)

xlabel('Time Window $T$','Interpreter','latex','FontSize',16); 
ylabel('$\| {\phi(t;\varphi)} \|^2_{L^2}$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',16)

filename = [pwd '/media/optimization/asymptotic_102_110_5K' ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])
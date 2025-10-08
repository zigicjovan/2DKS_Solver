%% RCG v RF eigenfunction validation

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

Jinitdata1 = load('Jinit_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_50.00.dat');
Joptdata1 = load('Jopt_s1_N_48_dt_0.005_K_1_1_L_1.02_1.10_T_0.10_50.00.dat');

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
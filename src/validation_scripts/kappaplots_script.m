        %{
        %%% timestep test %%%
        kappasinL = [ kappalist_48_1e2_sinL_p30 kappalist_48_5e3_sinL_p30 kappalist_48_1e3_sinL_p30  ] ;

        kappaerrorsinL= abs( 1 - kappasinL );

        % kappa test %
        h = figure;
        semilogx(logspace(-15,-1,15),kappasinL(:,1),'r--+')
        hold on
        semilogx(logspace(-15,-1,15),kappasinL(:,4),'r--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,7),'r--o')
        semilogx(logspace(-15,-1,15),kappasinL(:,2),'g--+')
        semilogx(logspace(-15,-1,15),kappasinL(:,5),'g--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,8),'g--o')
        semilogx(logspace(-15,-1,15),kappasinL(:,3),'b--+')
        semilogx(logspace(-15,-1,15),kappasinL(:,6),'b--*')
        semilogx(logspace(-15,-1,15),kappasinL(:,9),'b--o')

        yline(1,'--')
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylim([0.97 1.03])
        xlim([1e-15 1e-1])
        ylabel('$\kappa(\varepsilon)$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa test for $\varphi_{1}$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa1_file = [pwd '/data/kappa/kappasinL_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa1_file,'png');

        % kappa difference test %
        h = figure;
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,1),'r--+')
        hold on
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,4),'r--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,7),'r--o')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,2),'g--+')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,5),'g--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,8),'g--o')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,3),'b--+')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,6),'b--*')
        loglog(logspace(-15,-1,15),kappaerrorsinL(:,9),'b--o')

        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa difference test for $\varphi_{1}$", 'Interpreter','latex')
        legend("$N = 48, \Delta t = .01, \ell = 1.90", "$N = 48, \Delta t = .005, \ell = 1.90", "$N = 48, \Delta t = .001, \ell = 1.90", ...
            "$N = 48, \Delta t = .01, \ell = 2.25", "$N = 48, \Delta t = .005, \ell = 2.25", "$N = 48, \Delta t = .001, \ell = 2.25",  ...
            "$N = 48, \Delta t = .01, \ell = 2.60", "$N = 48, \Delta t = .005, \ell = 2.60", "$N = 48, \Delta t = .001, \ell = 2.60",  ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsinL_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa2_file,'png');
 
        %%% perturbation test %%%
        kappasinL30p = [ kappalist_48_1e2_sinL30_p1 kappalist_48_5e3_sinL30_p1 kappalist_48_1e2_sinL30_p10 kappalist_48_5e3_sinL30_p10  ] ;

        kappaerrorsinL30p = abs( 1 - kappasinL30p );

        i=0;
        % kappa difference test %
        i = i +2;
        h = figure;
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+1),'r--s')
        hold on
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+2),'r--*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+3),'r--o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+4),'m--s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+5),'m--*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+6),'m--o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+7),'g--s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+8),'g--*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+9),'g--o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+10),'b--s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+11),'b--*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i-1)*12+12),'b--o')

        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+1),'r:s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+2),'r:*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+3),'r:o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+4),'m:s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+5),'m:*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+6),'m:o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+7),'g:s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+8),'g:*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+9),'g:o')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+10),'b:s')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+11),'b:*')
        loglog(logspace(-15,-1,15),kappaerrorsinLp(:,(i)*12+12),'b:o')

        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
       
        xlim([1e-15 1e-1])
        %ylim([1e-4 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title("Kappa difference test for $\varphi_{1}, N=48, \varphi' = \varphi_{30}$", 'Interpreter','latex')
        legend("$\Delta t = .010, \ell = 1.15",  "$\Delta t = .010, \ell = 1.30", "$\Delta t = .010, \ell = 1.45", "$\Delta t = .010, \ell = 1.60", ...
            "$\Delta t = .010, \ell = 1.75", "$\Delta t = .010, \ell = 1.90", "$\Delta t = .010, \ell = 2.05", "$\Delta t = .010, \ell = 2.20", ...
            "$\Delta t = .010, \ell = 2.35", "$\Delta t = .010, \ell = 2.50",  "$\Delta t = .010, \ell = 2.65", "$\Delta t = .010, \ell = 2.80", ...
            "$\Delta t = .005, \ell = 1.15","$\Delta t = .005, \ell = 1.30",  "$\Delta t = .005, \ell = 1.45",  "$\Delta t = .005, \ell = 1.60", ...
            "$\Delta t = .005, \ell = 1.75",  "$\Delta t = .005, \ell = 1.90",  "$\Delta t = .005, \ell = 2.05",  "$\Delta t = .005, \ell = 2.20", ...
            "$\Delta t = .005, \ell = 2.35", "$\Delta t = .005, \ell = 2.50","$\Delta t = .005, \ell = 2.65", "$\Delta t = .005, \ell = 2.80", ...
            'Interpreter','latex', 'Location','southoutside', 'NumColumns',4)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsinLp30_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa2_file,'png');

        xlim([1e-10 1e-5])
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsinLp30zoom_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa2_file,'png');

        %%% time window test %%%

        kappaerr_T_1e2_sinL30_p30 = abs( 1 - kappalist_T_1e2_sinL30_p30 );

        i=0;
        % kappa difference test %
        i = i+1;
        h = figure;
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+1),'r--o')
        hold on
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+2),'r--x')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+3),'r--+')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+4),'r--*')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+5),'m--s')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+6),'m--d')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+7),'m--^')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+8),'m--p')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+9),'m--h')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+10),'g--<')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+11),'g--o')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+12),'g--x')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+13),'g--+')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+14),'g--*')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+15),'b--s')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+16),'b--d')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+17),'b--^')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+18),'b--p')
        loglog(logspace(-15,-1,15),kappaerr_T_1e2_sinL30_p30(:,(i-1)*length(timewindow)+19),'b--h')

        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        title('Kappa difference test for $\varphi_{30}: \ell = 2.80, N=48, \Delta t = .01, \varphi'' = \varphi_{30}$','Interpreter','latex');
        legend("$T=20", "$T=30$", "$T=40", "$T=50$", "$T=60", "$T=70$", "$T=80", "$T=90$", "$T=100$", ...
                "$T=110", "$T=120", "$T=130$", "$T=140", "$T=150$", "$T=160", "$T=170$", "$T=180", "$T=190$", "$T=200$", ...
            'Interpreter','latex', 'Location','eastoutside', 'NumColumns',1)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappaerrsinL30_ell280_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa2_file,'png');

        h = figure;
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((1-1)*length(timewindow)+1:(1)*length(timewindow)+0,1)),'r--o')
        hold on
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((2-1)*length(timewindow)+1:(2)*length(timewindow)+0,1)),'r--x')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((3-1)*length(timewindow)+1:(3)*length(timewindow)+0,1)),'r--+')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((4-1)*length(timewindow)+1:(4)*length(timewindow)+0,1)),'m--*')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((5-1)*length(timewindow)+1:(5)*length(timewindow)+0,1)),'m--s')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((6-1)*length(timewindow)+1:(6)*length(timewindow)+0,1)),'m--d')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((7-1)*length(timewindow)+1:(7)*length(timewindow)+0,1)),'g--o')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((8-1)*length(timewindow)+1:(8)*length(timewindow)+0,1)),'g--x')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((9-1)*length(timewindow)+1:(9)*length(timewindow)+0,1)),'g--+')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((10-1)*length(timewindow)+1:(10)*length(timewindow)+0,1)),'b--*')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((11-1)*length(timewindow)+1:(11)*length(timewindow)+0,1)),'b--s')
        semilogy(timewindow,abs(rieszlist_T_1e2_sinL30_p30((12-1)*length(timewindow)+1:(12)*length(timewindow)+0,1)),'b--d')
        xlabel('Time window $T$','Interpreter','latex'); 
        ylabel('$\langle \nabla\mathcal{J}_T (\varphi), \varphi''  \rangle_{L^2}$','Interpreter','latex');
        fontsize(12,"points")
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')    
        %title("Kappa test denominator for $\varphi_{s1}, N=48, \Delta t = .01, \varphi' = \varphi_{1}$", 'Interpreter','latex')
        title('Kappa test denominator','Interpreter','latex')
        subtitle(['$\varphi = \varphi_{s1}, \varphi'' = \varphi_{s1}, N = ' num2str(N) ', \Delta t = .01$'],'Interpreter','latex')        
        legend("$\ell = 1.15$", "$\ell = 1.30$", "$\ell = 1.45$", "$\ell = 1.60$", "$\ell = 1.75$", "$\ell = 1.90$", ...
            "$\ell = 2.05$", "$\ell = 2.20$", "$\ell = 2.35$", "$\ell = 2.50$", "$\ell = 2.65$", "$\ell = 2.80$", ...
            'Interpreter','latex', 'Location','northwest', 'NumColumns',3)
        frame = getframe(h);
        im = frame2im(frame);
        kappa2_file = [pwd '/data/kappa/kappadenomsinL_2DKS_N' num2str(N) '_dt' num2str(dt) '_T' num2str(T) '_lX' num2str(L_s1,'%.3f') '_lY' num2str(L_s2,'%.3f') '.png'];
        imwrite(im,kappa2_file,'png');
        %}


%{
kappalistnew = NaN(15,72);
rieszlistnew = NaN(72,1);
for i = 1:12
    kappalistnew(:,(0)*12+i) = kappalist2(:,(0)*12+i);
    kappalistnew(:,(1)*12+i) = kappalist1(:,(i-1)*3+1);
    kappalistnew(:,(2)*12+i) = kappalist1(:,(i-1)*3+2);
    kappalistnew(:,(3)*12+i) = kappalist2(:,(1)*12+i);
    kappalistnew(:,(4)*12+i) = kappalist1(:,(3)*12+(i-1)*3+1);
    kappalistnew(:,(5)*12+i) = kappalist1(:,(3)*12+(i-1)*3+2);
end
for i = 1:12
    rieszlistnew((0)*12+i,1) = rieszlist2((0)*12+i,1);
    rieszlistnew((1)*12+i,1) = rieszlist1((i-1)*3+1,1);
    rieszlistnew((2)*12+i,1) = rieszlist1((i-1)*3+2,1);
    rieszlistnew((3)*12+i,1) = rieszlist2((1)*12+i,1);
    rieszlistnew((4)*12+i,1) = rieszlist1((3)*12+(i-1)*3+1,1);
    rieszlistnew((5)*12+i,1) = rieszlist1((3)*12+(i-1)*3+2,1);
end
kappalist3 = [ kappalistnew(:,1:6), kappalistnew(:,10:18), kappalistnew(:,22:30) , kappalistnew(:,34:42) ...
 , kappalistnew(:,46:54), kappalistnew(:,58:66) , kappalistnew(:,70:72)] ;
rieszlist3 = [ rieszlistnew(1:6,1); rieszlistnew(10:18,1); rieszlistnew(22:30,1) ; rieszlistnew(34:42,1) ...
 ; rieszlistnew(46:54,1) ; rieszlistnew(58:66,1) ; rieszlistnew(70:72,1)] ;


plot(timestep,rieszlist3(1:3,1))
hold on
for i = 1:17
    plot(timestep,rieszlist3((i*3)+1:(i*3)+3,1))
end

legend('s1ps1 1.9 varying dt','s1ps1 2.25 varying dt','s1ps1 2.75 varying dt',...
    's1pstg1 1.9 varying dt','s1pstg1 2.25 varying dt','s1pstg1 2.75 varying dt',...
    's1pstg30 1.9 varying dt','s1pstg30 2.25 varying dt','s1pstg30 2.75 varying dt',...
    's30ps30 1.9 varying dt','s30ps30 2.25 varying dt','s30ps30 2.75 varying dt',...
    's30pstg1 1.9 varying dt','s30pstg1 2.25 varying dt','s30pstg1 2.75 varying dt',...
    's30pstg30 1.9 varying dt','s30pstg30 2.25 varying dt','s30pstg30 2.75 varying dt'...
    )

kappalist3 = abs(1-kappalist3);

%}


%%%%%%%%%%%%%%%%%%%%%

%kappalist = abs(1-kappalist);

kappalist = [ kappalist3(:,1:27), NaN(15,9) ];

j=1;
T=20;
% test 1: s1 s1
for i = 1:4
    figure
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+1),'r--x')
    hold on
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+2),'r--d')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+3),'r--o')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+4),'g--x')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+5),'g--d')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+6),'g--o')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+7),'b--x')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+8),'b--d')
    loglog(logspace(-15,-1,15),kappalist(:,(i-1)*9+9),'b--o')
    
    IC = strjoin(initialcondition(j),''); 
    pertIC = strjoin(kappapert(i),''); 
    N = 48;
    title('Kappa difference test','Interpreter','latex')
    subtitle(['$\varphi = \varphi_{' IC '}, \varphi'' = \varphi_{' pertIC '}, T = ' num2str(T,'%.0f') ', N = ' num2str(N) '$'],'Interpreter','latex')
    xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
    ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
    fontsize(12,"points")
    xlim([1e-15 1e-1])
    set(gca,'fontsize', 16) 
    set(gcf,'color','white')
    set(gca,'color','white')
    
    legend(['$\ell = ' num2str(L_scale(1),'%.3f') ', {\Delta}t = ' num2str(timestep(1),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(1),'%.3f') ', {\Delta}t = ' num2str(timestep(2),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(1),'%.3f') ', {\Delta}t = ' num2str(timestep(3),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(2),'%.3f') ', {\Delta}t = ' num2str(timestep(1),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(2),'%.3f') ', {\Delta}t = ' num2str(timestep(2),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(2),'%.3f') ', {\Delta}t = ' num2str(timestep(3),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(3),'%.3f') ', {\Delta}t = ' num2str(timestep(1),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(3),'%.3f') ', {\Delta}t = ' num2str(timestep(2),'%.4f') '$'],...
        ['$\ell = ' num2str(L_scale(3),'%.3f') ', {\Delta}t = ' num2str(timestep(3),'%.4f') '$'],...
        'Interpreter','latex','Location','north','NumColumns',3)
    hold off
end

%%%%%%%%%%%%

figure
semilogy(timewindow,rieszlist(1:2:19,1),'r-*')
hold on
semilogy(timewindow,rieszlist(21:2:39,1),'r-o')
semilogy(timewindow,rieszlist(2:2:20,1),'b-*')
semilogy(timewindow,rieszlist(22:2:40,1),'b-o')

dt = timestep(1);
N = 48;
title('Kappa test denominator','Interpreter','latex')
subtitle(['$\ell = ' num2str(L_scale(1),'%.3f') ', \varphi_{' strjoin(initialcondition(1),'') '}, \Delta t = ' num2str(dt,'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex')
xlabel('Time window $T$','Interpreter','latex'); 
ylabel('$\langle \nabla\mathcal{J}_T (\varphi), \varphi''  \rangle_{L^2}$','Interpreter','latex');
xlim([20 320])
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')

legend(['$\ell = ' num2str(L_scale(1),'%.3f') ', \varphi_{' strjoin(initialcondition(1),'') '}$'],...
    ['$\ell = ' num2str(L_scale(1),'%.3f') ',  \varphi_{' strjoin(initialcondition(2),'') '}$'],...
    ['$\ell = ' num2str(L_scale(2),'%.3f') ',  \varphi_{' strjoin(initialcondition(1),'') '}$'],...
    ['$\ell = ' num2str(L_scale(2),'%.3f') ',  \varphi_{' strjoin(initialcondition(2),'') '}$'],...
    'Interpreter','latex','Location','northwest')
hold off

kappalist = abs(1 - kappalist);

for j = 1:2
    for i = 1:2
        figure
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j),'r--x')
        hold on
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+2),'r--d')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+4),'r--s')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+6),'r--o')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+8),'g--x')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+10),'g--d')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+12),'g--o')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+14),'b--x')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+16),'b--d')
        loglog(logspace(-15,-1,15),kappalist(:,(i-1)*20+j+18),'b--o')
        
        IC = strjoin(initialcondition(j),''); 
        pertIC = IC; 
        N = 48;
        dt = timestep(1);
        title('Kappa difference test','Interpreter','latex')
        subtitle(['$\ell = ' num2str(L_scale(i)) ', \varphi = \varphi_{' IC '}, \varphi'' = \varphi_{' pertIC '}, {\Delta}t = ' num2str(timestep(1),'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex')
        xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
        ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
        fontsize(12,"points")
        xlim([1e-15 1e-1])
        set(gca,'fontsize', 16) 
        set(gcf,'color','white')
        set(gca,'color','white')
        
        legend(['$T = ' num2str(timewindow(1),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(2),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(3),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(4),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(5),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(6),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(7),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(8),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(9),'%.0f') '$'],...
                ['$T = ' num2str(timewindow(10),'%.0f') '$'],...
            'Interpreter','latex','Location','eastoutside','NumColumns',1)
        hold off
    end
end

%%%%%%%%%%%%

semilogy(20:20:500,rieszlists1236,'r-*')
title('Kappa test denominator','Interpreter','latex')
xlim([20 320])
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')
xlabel('Time window $T$','Interpreter','latex'); 
ylabel('$\langle \nabla\mathcal{J}_T (\varphi), \varphi''  \rangle_{L^2}$','Interpreter','latex');
subtitle(['$\ell = ' num2str(L_scale(1),'%.3f') ', \varphi_{' strjoin(initialcondition(1),'') '}, \Delta t = ' num2str(dt,'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex')
subtitle(['$\ell = ' num2str(L_scale(1),'%.3f') ', \varphi_{' strjoin(initialcondition(1),'') '}, \Delta t = ' num2str(dt,'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex')

figure
%epscheck = logspace(-15,-1,15);
epscheck = [.001,.0025,.005,.0075,.01,.025,.05,.075,.1,.25,.5,.75,1];
loglog(epscheck,kappalist(:,11),'r--x')
hold on
loglog(epscheck,kappalist(:,12),'r--d')
loglog(epscheck,kappalist(:,13),'r--s')
loglog(epscheck,kappalist(:,14),'r--o')
loglog(epscheck,kappalist(:,15),'g--x')
loglog(epscheck,kappalist(:,16),'g--d')
loglog(epscheck,kappalist(:,17),'g--s')
loglog(epscheck,kappalist(:,18),'b--x')
loglog(epscheck,kappalist(:,19),'b--d')
loglog(epscheck,kappalist(:,20),'b--s')
%loglog(epscheck,kappalists1236(:,11),'b--x')
%loglog(epscheck,kappalists1236(:,12),'b--o')
%loglog(epscheck,kappalists1236(:,13),'b--s')
%loglog(epscheck,kappalists1236(:,14),'b--d')
%loglog(epscheck,kappalists1236(:,15),'b--+')
%loglog(epscheck,kappalists1236(:,16),'k--x')

IC = strjoin(initialcondition(1),''); 
pertIC = IC; 
N = 48;
dt = timestep(1);
title('Kappa difference test','Interpreter','latex')
subtitle(['$\ell = ' num2str(L_scale(2)) ', \varphi = \varphi_{' IC '}, \varphi'' = \varphi_{' pertIC '}, {\Delta}t = ' num2str(timestep(1),'%.4f') ', N = ' num2str(N) '$'],'Interpreter','latex')
xlabel('Perturbation magnitude $\varepsilon$','Interpreter','latex'); 
ylabel('$| 1 - \kappa(\varepsilon)|$','Interpreter','latex');
fontsize(12,"points")
%xlim([1e-15 1e-1])
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')

legend(['$T = ' num2str(timewindow(1),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(2),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(3),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(4),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(5),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(6),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(7),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(8),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(9),'%.0f') '$'],...
        ['$T = ' num2str(timewindow(10),'%.0f') '$'],...
    'Interpreter','latex','Location','eastoutside','NumColumns',1)
hold off

% max kappa plot
h = figure;
set(gcf,'Position',[100 100 900 750])
plot(L_scale,maxkappalist(:,1),'r-*')
hold on
plot(L_scale,maxkappalist(:,2),'b-o')
ylim([0 400])
xline(1.48,'--',{'1.48'},'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
xlim([L_scale(1) L_scale(end)])
xlabel('Domain factor $\ell$','Interpreter','latex'); 
ylabel('Time window $T$','Interpreter','latex');
%fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title('Maximum feasible optimization time windows','Interpreter','latex')
subtitle(['$\varphi = \varphi_{s1}, {\Delta}t = ' num2str(0.005) ', N = ' num2str(48) '$'],'Interpreter','latex','FontSize',14)
legend( '$\varphi'' = \varphi_{s1}$', '$\varphi'' = \varphi_{stg30}$' ,'Interpreter','latex','Location','northeast')
filename = [pwd '/media/kappa/maxkappa'];
exportgraphics(h,[filename '.pdf'])
saveas(h,[filename '.fig'])

% max riesz plot
h = figure;
set(gcf,'Position',[100 100 900 750])
semilogy(L_scale,maxrieszlist(:,1),'r-*')
hold on
semilogy(L_scale,maxrieszlist(:,2),'b-o')
%ylim([0 400])
xline(1.48,'--',{'1.48'},'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
xlim([L_scale(1) L_scale(end)])
xlabel('Domain factor $\ell$','Interpreter','latex'); 
ylabel('$\langle \nabla\mathcal{J}_T (\varphi), \varphi''  \rangle_{L^2}$','Interpreter','latex');
title('Magnitude of $\kappa(\varepsilon)$ denominator for maximum feasible time windows','Interpreter','latex')
subtitle(['$\varphi = \varphi_{s1}, {\Delta}t = ' num2str(0.005) ', N = ' num2str(48) '$'],'Interpreter','latex','FontSize',14)
%fontsize(12,"points")
legend( '$\varphi'' = \varphi_{s1}$', '$\varphi'' = \varphi_{stg30}$' ,'Interpreter','latex','Location','northwest')
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
filename = [pwd '/media/kappa/maxriesz'];
exportgraphics(h,[filename '.pdf'])
saveas(h,[filename '.fig'])

% 2.36 and 2.78 trajectories
h = figure;
set(gcf,'Position',[100 100 900 750])
semilogy(timewindow,s136,'r-')
hold on
semilogy(timewindow,stg136,'b-')
%ylim([0 400])
%xline(1.48,'--',{'1.48'},'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
xlim([0 T])
xlabel('Time $t$','Interpreter','latex'); 
ylabel('$||{\phi(t;\varphi)}||_{L^2}$','Interpreter','latex');
title('Evolution of $L^2$ norm','Interpreter','latex')
subtitle(['$\ell = 2.36$, ${\Delta}t = ' num2str(0.005) ', N = ' num2str(48) '$'],'Interpreter','latex','FontSize',14)
%fontsize(12,"points")
legend( '$\varphi = \varphi_{s1}$', '$\varphi = \varphi_{stg1}$' ,'Interpreter','latex','Location','southeast')
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
filename = [pwd '/media/energy/s1stg1_236.pdf'];
exportgraphics(h,filename)

h = figure;
set(gcf,'Position',[100 100 900 750])
semilogy(timewindow,s178,'r-')
hold on
semilogy(timewindow,stg178,'b-')
%ylim([0 400])
%xline(1.48,'--',{'1.48'},'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
xlim([0 T])
xlabel('Time $t$','Interpreter','latex'); 
ylabel('$||{\phi(t;\varphi)}||_{L^2}$','Interpreter','latex');
title('Evolution of $L^2$ norm','Interpreter','latex')
subtitle(['$\ell = 2.78$, ${\Delta}t = ' num2str(0.005) ', N = ' num2str(48) '$'],'Interpreter','latex','FontSize',14)
%fontsize(12,"points")
legend( '$\varphi = \varphi_{s1}$', '$\varphi = \varphi_{stg1}$' ,'Interpreter','latex','Location','southeast')
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
filename = [pwd '/media/energy/s1stg1_278.pdf'];
exportgraphics(h,filename)
function maxL2inT = plot_2DKS(save_each, solplot, IC, N, dt, T, K, L_s1, L_s2, Ntime_save_max, utility1, utility2)

if not(isfolder([pwd  '/data/energyL2' ]))                                           % create local directories for data storage
    mkdir([pwd  '/data/energyL2' ])
    mkdir([pwd  '/data/energyL2_t' ]);
    mkdir([pwd  '/data/spectrum' ]);
    mkdir([pwd  '/data/kappa' ]);
    mkdir([pwd  '/media' ]);
    mkdir([pwd  '/media/movies' ]);
    mkdir([pwd  '/media/energy' ]);
    mkdir([pwd  '/media/optimization' ]);
    mkdir([pwd  '/media/figures' ]);
    mkdir([pwd  '/media/figures/state' ]);
    mkdir([pwd  '/media/kappa' ]);
    addpath([pwd  '/data/energyL2' ])
    addpath([pwd  '/data/energyL2_t' ]);
    addpath([pwd  '/data/spectrum' ]);
    addpath([pwd  '/data/kappa' ]);
    addpath([pwd  '/media' ]);
    addpath([pwd  '/media/movies' ]);
    addpath([pwd  '/media/energy' ]);
    addpath([pwd  '/media/optimization' ]);
    addpath([pwd  '/media/figures' ]);
    addpath([pwd  '/media/figures/state' ]);
    addpath([pwd  '/media/kappa' ]);
end

parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
parfiglist = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];
optparfiglist = ['$\varphi = \tilde{\varphi}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];
originalIC = utility1;
tol = utility2(1);
if length(utility2) > 1
    if utility2(2) == 1
        optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) '_RCG' ];
    elseif utility2(2) == 0
        optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) '_RG' ];
    end
else
    optparameters = [ originalIC '_' parameterlist '_tol_' num2str(tol) ];
end

% number of timesteps
time1 = ceil(T/dt); 
time2 = ceil(time1/save_each);
if T >= 0
    Ntime = max(time1,time2);
    Ntime_save = min(time1,time2);
    save_each = ceil(Ntime/Ntime_save);
    Ntime = ceil(Ntime/save_each);
end

switch solplot
    case {'norms','gif','diagnostics','initial','terminal','optdiag'}
        savedata = 1;
        [maxL2inT,u_IC,u_TC,energyL2,energyH1,energyH2,astripwidth,v_mean,projcoeffradialevolution,projcoeffmodeevolution] = ...
            process_energy(savedata,IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
            parameterlist,optparameters,parfiglist,optparfiglist);
end

switch IC
    case {'optimized'}
        savedata = 0;
        [~,u_IC_og,u_TC_og,energyL2_og,energyH1_og,energyH2_og,astripwidth_og,v_mean_og,projcoeffradialevolution_og,projcoeffmodeevolution_og] = ...
            process_energy(savedata,originalIC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
            parameterlist,optparameters,parfiglist,optparfiglist);
end

switch solplot
    case 'gif'

        process_gif(IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
                    energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                    parameterlist,optparameters,parfiglist,optparfiglist);

    case 'diagnostics'

        process_figure('diagnostics',originalIC, IC, dt, T, N, K, L_s1, L_s2, energyL2_og,utility2,Ntime,Ntime_save_max,tol,... 
                energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                parameterlist,optparameters,parfiglist,optparfiglist)

    case 'state'

        process_figure('state',originalIC, IC, dt, T, N, K, L_s1, L_s2, u_IC,u_TC,Ntime,Ntime_save_max,tol,... 
                energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                parameterlist,optparameters,parfiglist,optparfiglist)

    case 'kappa'

        kappa = utility1;
        pertIC = utility2;
        kappaerror = abs( 1 - kappa );

        h = figure('Visible', 'off');
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
        subtitle(['$\varphi'' = \varphi_{' pertIC '}$, ' parfiglist],'Interpreter','latex','FontSize',14)
        filename = [pwd '/media/kappa/kappaerr_p' pertIC '_' parameterlist '.pdf'];
        exportgraphics(h,filename)

        kappa_file = [pwd '/data/kappa/kappa_p' pertIC '_' parameterlist '.dat'];
        writematrix(kappaerror, kappa_file,'Delimiter','tab');
    case 'optdiag'
        
        close all
        %% opt gif

        process_gif(IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
                    energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                    parameterlist,optparameters,parfiglist,optparfiglist);
        fprintf('Saved optimized evolution video at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))

        close all
        if utility2(4) > 0
            %% diagnostics

            process_figure('energy',originalIC, IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,tol,... 
                energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                parameterlist,optparameters,parfiglist,optparfiglist)

            process_figure('diagnostics',originalIC, IC, dt, T, N, K, L_s1, L_s2, energyL2_og,utility2,Ntime,Ntime_save_max,tol,... 
                energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                parameterlist,optparameters,parfiglist,optparfiglist)
        end

	    fprintf('Saved optimization diagnostic figures at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))

        %% opt states
        process_figure('state',originalIC, IC, dt, T, N, K, L_s1, L_s2, u_IC,u_TC,Ntime,Ntime_save_max,tol,... 
                energyL2,energyH1,energyH2,v_mean,astripwidth,projcoeffradialevolution,projcoeffmodeevolution,...
                parameterlist,optparameters,parfiglist,optparfiglist)
        
        %% finish opt plots
        IC = originalIC;
        parameterlist = [IC '_N_' num2str(N) '_dt_' num2str(dt) '_K_' num2str(K,'%.0f') '_Ls1_' num2str(L_s1,'%.2f') '_Ls2_' num2str(L_s2,'%.2f') '_T_' num2str(T) ];
        parfiglist = ['$\varphi = \varphi_{' IC '}, N = ' num2str(N) ', {\Delta}t = ' num2str(dt) ', K = ' num2str(K,'%.0f') ', L_1 = 2\pi(' num2str(L_s1,'%.2f') '), L_2 = 2\pi(' num2str(L_s2,'%.2f') '), T = ' num2str(T,'%.2f') '$'];

        %% original states
        process_figure('state',originalIC, IC, dt, T, N, K, L_s1, L_s2, u_IC_og,u_TC_og,Ntime,Ntime_save_max,tol,... 
                energyL2_og,energyH1_og,energyH2_og,v_mean_og,astripwidth_og,projcoeffradialevolution_og,projcoeffmodeevolution_og,...
                parameterlist,optparameters,parfiglist,optparfiglist)

        %% original gif
        
        process_gif(IC, dt, T, N, K, L_s1, L_s2, utility1,utility2,Ntime,Ntime_save_max,... 
                    energyL2_og,energyH1_og,energyH2_og,v_mean_og,astripwidth_og,projcoeffradialevolution_og,projcoeffmodeevolution_og,...
                    parameterlist,optparameters,parfiglist,optparfiglist);
        fprintf('Saved original evolution video at %01dh%02dm%02ds\n',floor(toc/3600),floor(mod(toc/60,60)),floor(mod(toc,60)))

end

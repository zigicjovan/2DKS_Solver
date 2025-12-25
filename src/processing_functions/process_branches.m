function process_branches(specifiedrange,Kstart,Kend,Knum,ellstart,ellend,ellgap)                                    

%% identify potential range (specifiedrange indicates whether to use it or plot all)
Krange = round(logspace(Kstart,Kend,Knum));                             
ellrange = ellstart:ellgap:ellend;        

%% sort all files by T and collect argmax Jopt(T)
branchDir = [pwd '/data/branches'];
branchfiles = dir(fullfile(branchDir, 'Jopt_*.dat'));

maxJopt = NaN(numel(branchfiles),6);
branch_groups = cell(numel(branchfiles),3);

for k = 1:numel(branchfiles)
    fname = fullfile(branchfiles(k).folder, branchfiles(k).name);
    branch = readmatrix(fname);
    branch = sortrows(branch, 3);
    [~,idx] = max(branch(:,4));
    maxJopt(k,:) = branch(idx,:);
    branch_groups{k,1} = branch(1,1);
    branch_groups{k,2} = branch(1,2);
    branch_groups{k,3} = branch;
    writematrix(branch, fname, 'Delimiter', 'tab');
end

maxJopt = sortrows(maxJopt, [2,1]);
branch_groups = sortrows(branch_groups, [2,1]);

%% plot entire or specified range 
if specifiedrange == 1
    maxJopt(maxJopt(:,1) > 10^Kend, :) = [];
    maxJopt(maxJopt(:,2) > ellend, :) = [];
    branch_groups(cell2mat(branch_groups(:,1)) > 10^Kend, :) = [];
    branch_groups(cell2mat(branch_groups(:,2)) > ellend, :) = [];
end

%% group ell by max J(T) and argmax J(T)
tol = 1e-10;
ell_vals = unique(round(maxJopt(:,2)/tol)*tol);
K_vals = unique(round(maxJopt(:,1)/tol)*tol);

ell_groups = cell(numel(ell_vals),2);
K_groups = [ [0 , ell_vals'] ; [ K_vals , NaN(numel(K_vals),numel(ell_vals)) ] ];
for i = 1:numel(ell_vals)
    ell_groups{i,1} = ell_vals(i);
    ell_groups{i,2} = maxJopt(abs(maxJopt(:,2) - ell_vals(i)) < tol, :);
    T_vals = maxJopt(abs(maxJopt(:,2) - ell_vals(i)) < tol, 3);
    K_groups(2:(size(T_vals,1) + 1),i+1) = T_vals;
end

%% argmax J(T) figure with exponential law

h = figure;
T_Legend = [];     % handles to include
T_legendlabels  = {};     % legend labels
T_loglinfits = NaN(numel(K_vals),2);

for k = 1:numel(K_vals)
    T_fitlength = length(rmmissing(K_groups(k+1,2:end)));
    T_loglinfits(k,:) = polyfit(K_groups(1,2:T_fitlength),log10(K_groups(k+1,2:T_fitlength)), 1);
    curveplot = semilogy(K_groups(1,2:end),(K_groups(k+1,2:end)),'-o','MarkerSize',3);
    T_Legend(end+1) = curveplot;
    T_legendlabels{end+1}  = ['$' num2str(K_groups(k+1,1),'%.0f') ' \approx 10^{' num2str(T_loglinfits(k,1),'%.3f') '\ell ' num2str(T_loglinfits(k,2),'%.3f') '}$'];
    hold on
end
T_loglinfits = rmmissing(T_loglinfits);
avgSlope = mean(T_loglinfits(:,1));
stdSlope = std(T_loglinfits(:,1));
coeffmin = T_loglinfits(end,2);
coeffmax = T_loglinfits(1,2);
powerlawmin = min(T_loglinfits(:,1));
powerlawmax = max(T_loglinfits(:,1));
curveplot = semilogy(ellrange, (10.^(coeffmin + ellrange.*powerlawmin)),'b--','MarkerSize',3);
T_Legend(end+1) = curveplot;
T_legendlabels{end+1}  = ['$10^{' num2str(powerlawmin,'%.3f') '\ell ' num2str(coeffmin,'%.3f') '}$'];
curveplot = semilogy(ellrange, (10.^(coeffmax + ellrange.*powerlawmax)),'r--','MarkerSize',3);
T_Legend(end+1) = curveplot;
T_legendlabels{end+1}  = ['$10^{' num2str(powerlawmax,'%.3f') '\ell ' num2str(coeffmax,'%.3f') '}$'];

hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1250 700])
xlabel('Domain size $\ell$','Interpreter','latex'); 
ylabel('Peak Time Window $\tilde{T}_L$','Interpreter','latex');
xlim([ellstart ellend])
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = ['$K = \left[10^{' num2str(Kstart) '},10^{' num2str(Kstart+0.5) '}, 10^{' num2str(Kend) '}\right], \ell = [' num2str(ellstart) ',' num2str(ellstart+ellgap) ', ... , ' num2str(ellend) ']' ...
    ',  \tilde{T}_{L} \approx 10^{(' num2str(avgSlope,'%.3f') ' \pm ' num2str(stdSlope,'%.3f') ')\ell}$'];
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
nCols = ceil(numel(T_Legend) / 5);
legend(T_Legend, T_legendlabels,'Interpreter','latex','Location','southeast','NumColumns',1,'Box','off','FontSize',12)

filename = [pwd '/media/optimization/argmaxJ_T_K_' num2str(Kstart) '_' num2str(Kend) '_L_' num2str(ellstart) '_' num2str(ellend)  ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% branches figure
branchcolor = 0;
h = figure;
K_Legend = [];     % handles to include
K_legendlabels  = {};     % legend labels

for i = 1:size(branch_groups,1)
    if i > 1 && branch_groups{i,2} ~= branch_groups{i-1,2}
        branchcolor = branchcolor + 1;
    end
    currentbranch = branch_groups{i,3};
    if mod(branchcolor,5) == 0
        curveplot = loglog(currentbranch(:,3),currentbranch(:,4),'b-o','MarkerSize',3);
    elseif mod(branchcolor,5) == 1
        curveplot = loglog(currentbranch(:,3),currentbranch(:,4),'r-o','MarkerSize',3);
    elseif mod(branchcolor,5) == 2
        curveplot = loglog(currentbranch(:,3),currentbranch(:,4),'g-o','MarkerSize',3);
    elseif mod(branchcolor,5) == 3
        curveplot = loglog(currentbranch(:,3),currentbranch(:,4),'m-o','MarkerSize',3);
    elseif mod(branchcolor,5) == 4
        curveplot = loglog(currentbranch(:,3),currentbranch(:,4),'k-o','MarkerSize',3);
    end
    hold on
    if (i == 1) || (i > 1 && branch_groups{i,2} ~= branch_groups{i-1,2})
        K_Legend(end+1) = curveplot;
        K_legendlabels{end+1}  = sprintf('%.2f', branch_groups{i,2});
    end
end
hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1250 700])
xlabel('Time Window $T$','Interpreter','latex'); 
ylabel('$\| {\phi(T;\varphi)} \|^2_{L^2}$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = ['$K = \left[10^{' num2str(Kstart) '},10^{' num2str(Kstart+0.5) '}, 10^{' num2str(Kend) '}\right], \ell = [' num2str(ellstart) ',' num2str(ellstart+ellgap) ', ... , ' num2str(ellend) ']$'];
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
nCols = ceil(numel(K_Legend) / 5);
legend(K_Legend, K_legendlabels,'Interpreter','latex', 'Location', 'southwest','NumColumns',nCols,'Box','off');

filename = [pwd '/media/optimization/branches_K_' num2str(Kstart) '_' num2str(Kend) '_L_' num2str(ellstart) '_' num2str(ellend)  ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% max J(T) figure with power law
h = figure;
K_Legend = [];     % handles to include
K_legendlabels  = {};     % legend labels
K_loglogfits = NaN(size(ell_groups,1),2); % for LLS fit

for i = 1:size(ell_groups,1)
    currentell = ell_groups{i,2};
    curveplot = loglog(currentell(:,1),currentell(:,4),'-o','MarkerSize',3);
    if size(currentell,1) > 1
        loglogfit_i = polyfit(log10(currentell(:,1)),log10(currentell(:,4)), 1);
        K_loglogfits(i,:) = loglogfit_i;
    end
    K_Legend(end+1) = curveplot;
    K_legendlabels{end+1}  = ['$' num2str(ell_groups{i,1},'%.2f') ' \approx 10^{' num2str(loglogfit_i(2),'%.3f') '}K^{' num2str(loglogfit_i(1),'%.3f') '}$'];
    hold on
end
K_loglogfits = rmmissing(K_loglogfits);
avgSlope = mean(K_loglogfits(:,1));
stdSlope = std(K_loglogfits(:,1));
coeffmin = K_loglogfits(1,2);
coeffmax = K_loglogfits(end,2);
powerlawmin = min(K_loglogfits(:,1));
powerlawmax = max(K_loglogfits(:,1));
curveplot = loglog(Krange, 10^(coeffmin).*(Krange.^(powerlawmin)),'b--','MarkerSize',3);
K_Legend(end+1) = curveplot;
K_legendlabels{end+1}  = ['$10^{' num2str(coeffmin,'%.3f') '} K^{' num2str(powerlawmin,'%.3f') '}$'];
curveplot = loglog(Krange, 10^(coeffmax).*(Krange.^(powerlawmax)),'r--','MarkerSize',3);
K_Legend(end+1) = curveplot;
K_legendlabels{end+1}  = ['$10^{' num2str(coeffmax,'%.3f') '} K^{' num2str(powerlawmax,'%.3f') '}$'];

hold off

set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white') 
set(gcf,'Position',[100 100 1250 700])
xlabel('Initial Energy $K$','Interpreter','latex'); 
ylabel('Peak Terminal Energy $\tilde{K}_L$','Interpreter','latex');
title('Optimized Finite-Time $L^2$ energy for 2D Kuramoto-Sivashinsky','Interpreter','latex')
parfiglistInterval = ['$K = \left[10^{' num2str(Kstart) '},10^{' num2str(Kstart+0.5) '}, 10^{' num2str(Kend) '}\right], \ell = [' num2str(ellstart) ',' num2str(ellstart+ellgap) ', ... , ' num2str(ellend) ']' ...
    ',  \tilde{K}_L \approx K^{' num2str(avgSlope,'%.3f') ' \pm ' num2str(stdSlope,'%.3f') '}$'];
subtitle(parfiglistInterval,'Interpreter','latex','FontSize',14)
nCols = ceil(numel(K_Legend) / 5);
legend(K_Legend, K_legendlabels,'Interpreter','latex','Location','southoutside','NumColumns',5,'Box','off','FontSize',12)

filename = [pwd '/media/optimization/maxJopt_K_' num2str(Kstart) '_' num2str(Kend) '_L_' num2str(ellstart) '_' num2str(ellend)  ];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])


%% max T law


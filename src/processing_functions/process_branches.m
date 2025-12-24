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

%% plot all (or specified) metrics

if specifiedrange == 0

    %% group by ell
    tol = 1e-10;
    ell_vals = unique(round(maxJopt(:,2)/tol)*tol);
    
    ell_groups = cell(numel(ell_vals),2);
    for i = 1:numel(ell_vals)
        ell_groups{i,1} = ell_vals(i);
        ell_groups{i,2} = maxJopt(abs(maxJopt(:,2) - ell_vals(i)) < tol, :);
    end

    %% branches figure
    branchcolor = 0;
    h = figure;
    hLegend = [];     % handles to include
    legendlabels  = {};     % legend labels

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
            hLegend(end+1) = curveplot;
            legendlabels{end+1}  = sprintf('%.2f', branch_groups{i,2});
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
    nCols = ceil(numel(hLegend) / 5);
    legend(hLegend, legendlabels, 'Location', 'southwest','NumColumns',nCols,'Box','off');

    filename = [pwd '/media/optimization/branches_K_' num2str(Kstart) '_' num2str(Kend) '_L_' num2str(ellstart) '_' num2str(ellend)  ];
    saveas(h,[filename '.fig'])
    exportgraphics(h,[filename '.pdf'])

    %% maxJopt figure
    h = figure;
    hLegend = [];     % handles to include
    legendlabels  = {};     % legend labels
    linearfits = NaN(size(ell_groups,1),2); % for LLS fit

    for i = 1:size(ell_groups,1)
        currentell = ell_groups{i,2};
        curveplot = loglog(currentell(:,1),currentell(:,4),'-o','MarkerSize',3);
        if size(currentell,1) > 1
            linearfit_i = polyfit(log10(currentell(:,1)),log10(currentell(:,4)), 1);
            linearfits(i,:) = linearfit_i;
        end
        hLegend(end+1) = curveplot;
        legendlabels{end+1}  = ['$' num2str(ell_groups{i,1},'%.3f') ' \approx 10^{' num2str(linearfit_i(2),'%.3f') '}K^{' num2str(linearfit_i(1),'%.3f') '}$'];
        hold on
    end
    linearfits = rmmissing(linearfits);
    avgSlope = mean(linearfits(:,1));
    stdSlope = std(linearfits(:,1));
    coeffmin = linearfits(1,2);
    coeffmax = linearfits(end,2);
    xlin = logspace(Kstart,Kend,Knum);
    powerlawmin = min(linearfits(:,1));
    powerlawmax = max(linearfits(:,1));
    curveplot = loglog(xlin, 10^(coeffmin).*(xlin.^(powerlawmin)),'b--','MarkerSize',3);
    hLegend(end+1) = curveplot;
    legendlabels{end+1}  = ['$10^{' num2str(coeffmin,'%.3f') '} K^{' num2str(powerlawmin,'%.3f') '}$'];
    curveplot = loglog(xlin, 10^(coeffmax).*(xlin.^(powerlawmax)),'r--','MarkerSize',3);
    hLegend(end+1) = curveplot;
    legendlabels{end+1}  = ['$10^{' num2str(coeffmax,'%.3f') '} K^{' num2str(powerlawmax,'%.3f') '}$'];

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
    nCols = ceil(numel(hLegend) / 5);
    legend(hLegend, legendlabels,'Interpreter','latex','Location','southoutside','NumColumns',5,'Box','off','FontSize',12)

    filename = [pwd '/media/optimization/maxJopt_K_' num2str(Kstart) '_' num2str(Kend) '_L_' num2str(ellstart) '_' num2str(ellend)  ];
    saveas(h,[filename '.fig'])
    exportgraphics(h,[filename '.pdf'])
else
    
end

%% max T law
numell = size(ell_groups,1);                    % number of rows
maxTlaw = NaN(Knum,numell);                     % one vector per row

for i = 1:numell
    maxK = ell_groups{i,2};                     % numeric matrix in column 2
    numK = size(maxK,1);
    maxTlaw(1:numK,i) = maxK(:,3);              % extract T value of max J
end

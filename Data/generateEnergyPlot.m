clear;
close all;
clc;

% Run this script from the directory containing the .dat files.

%% Common parameters

ells = [1.05 1.20 1.35 1.50 1.65 1.80 1.95];
K_values = [1.0e4 3.2e4 1.0e5 3.2e5 1.0e6];
K_file_strings = { ...
    '1.0e+04', ...
    '3.2e+04', ...
    '1.0e+05', ...
    '3.2e+05', ...
    '1.0e+06'};
K_labels = { ...
    '10^{4.0}', ...
    '10^{4.5}', ...
    '10^{5.0}', ...
    '10^{5.5}', ...
    '10^{6.0}'};

fontSize = 20;
figPosition = [100 100 1200 800];

% MATLAB default color order gives a consistent distinct color per curve.
colors = lines(max(numel(ells), numel(K_values)));

%% ========================================================================
%  ENERGY EVOLUTION
% ========================================================================

%% Create figures first
figL2 = figure('Position', figPosition);
axL2 = axes(figL2);
hold(axL2, 'on');

figH1 = figure('Position', figPosition);
axH1 = axes(figH1);
hold(axH1, 'on');

figH2 = figure('Position', figPosition);
axH2 = axes(figH2);
hold(axH2, 'on');

for j = numel(K_values):-1:1
    for i = numel(ells):-1:1
    
        ell = ells(i);
        Kstr = K_file_strings{j};
    
        filename = sprintf('*K_%s_ell1_%.2f_ell2_%.2f_*.dat', Kstr, ell, ell);
    
        try
            fname = dir(filename).name;
            fid = fopen(fname, 'r');
            
            % Pass 1: count text rows
            lengthdata = 0;
            while ischar(fgetl(fid))
                lengthdata = lengthdata + 1;
            end
            samplestep = max(1, round(lengthdata / 1e4));
            
            % Pass 2: read only every samplestep-th row
            frewind(fid);  
            nSamples = ceil(lengthdata / samplestep);
            timept = zeros(nSamples, 1);
            energyL2 = zeros(nSamples, 1);
            energyH1 = zeros(nSamples, 1);
            energyH2 = zeros(nSamples, 1);
            
            ik = 1;
            jk = 1;
            
            while true
                line = fgetl(fid);
                if ~ischar(line)
                    break
                end
            
                if mod(ik-1, samplestep) == 0
                    vals = sscanf(line, '%f');
                    timept(jk) = vals(1);
                    energyL2(jk) = vals(2);
                    energyH1(jk) = vals(3);
                    energyH2(jk) = vals(4);
                    jk = jk + 1;
                end
            
                ik = ik + 1;
            end
            
            fclose(fid);
            
            %% Remove unused preallocated entries
            timept   = timept(1:jk-1);
            energyL2 = energyL2(1:jk-1);
            energyH1 = energyH1(1:jk-1);
            energyH2 = energyH2(1:jk-1);

            %% Plot all three
            loglog(axL2, timept, energyL2, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
            loglog(axH1, timept, energyH1, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
            loglog(axH2, timept, energyH2, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        catch
            if exist('fid','var') && fid ~= -1
                fclose(fid);
            end
        end
    end
end

%% Common formatting
axesList = [axL2, axH1, axH2];
for ax = axesList
    grid(ax, 'on');
    box(ax, 'on');
    %xlim(ax, 'tight');
    xlim(ax, [1e-6 1.2e0])
    ylim(ax, 'tight');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
    xlabel(ax, 'Time $t$', 'Interpreter', 'latex', 'FontSize', fontSize);
end

%% L2
ylabel(axL2, '$L^2$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{L^2}$', ...
    'Interpreter', 'latex', 'FontSize', fontSize);
title(axL2, '$L^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', 'FontSize', fontSize);

%% H1
ylabel(axH1, '$H^1$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{H^1}$', ...
    'Interpreter', 'latex', 'FontSize', fontSize);
title(axH1, '$H^1$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', 'FontSize', fontSize);

%% H2
ylabel(axH2, '$H^2$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{H^2}$', ...
    'Interpreter', 'latex', 'FontSize', fontSize);
title(axH2, '$H^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', 'FontSize', fontSize);

%% Same subtitle
subtitleText = '$K=[10^4,10^{4.5},\ldots,10^6]$ and $\ell=[1.05,1.20,\ldots,1.95]$';

subtitle(axL2, subtitleText, 'Interpreter','latex', 'FontSize',fontSize);
subtitle(axH1, subtitleText, 'Interpreter','latex', 'FontSize',fontSize);
subtitle(axH2, subtitleText, 'Interpreter','latex', 'FontSize',fontSize);

%% Export
exportgraphics(figL2, 'energyL2evolution.pdf', 'ContentType', 'vector');
exportgraphics(figH1, 'energyH1evolution.pdf', 'ContentType', 'vector');
exportgraphics(figH2, 'energyH2evolution.pdf', 'ContentType', 'vector');

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
%  1. L^2 ENERGY EVOLUTION
% ========================================================================

figure('Position', figPosition);
hold on;

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
            energy = zeros(nSamples, 1);
            
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
                    energy(jk) = vals(2);
                    jk = jk + 1;
                end
            
                ik = ik + 1;
            end
            
            fclose(fid);
            
            timept = timept(1:jk-1);
            energy = energy(1:jk-1);
            loglog( timept, energy, '-', ...
                    'Color', colors(i,:), ...
                    'LineWidth', 1.5, ...
                    'MarkerSize', 6);
        catch
        end
    end
end
xlabel('Time $t$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

ylabel('$L^2$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{L^2}$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

title('$L^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
      'Interpreter', 'latex', ...
      'FontSize', fontSize);

subtitle([ '$K=[10^4,10^{4.5},\ldots,10^6]$ and ' ...
    '$\ell=[1.05,1.20,\ldots,1.95]$'], ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;
%xlim([1e4 1e6]);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

exportgraphics(gcf, 'energyL2evolution.pdf', 'ContentType', 'vector');

%% ========================================================================
%  2. H^1 ENERGY EVOLUTION
% ========================================================================

figure('Position', figPosition);
hold on;

for j = numel(K_values):-1:1
    for i = numel(ells):-1:1
    
        ell = ells(i);
        Kstr = K_file_strings{j};
    
        filename = sprintf('*K_%s_ell1_%.2f_ell2_%.2f_*.dat', Kstr, ell, ell);
    
        try
            data = readmatrix(dir(filename).name);
            timept = data(:,1);
            energy = data(:,3);
            loglog( timept, energy, '-', ...
                    'Color', colors(i,:), ...
                    'LineWidth', 1.5, ...
                    'MarkerSize', 6);
        catch
        end
    end
end
xlabel('Time $t$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

ylabel('$H^1$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{H^1}$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

title('$H^1$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
      'Interpreter', 'latex', ...
      'FontSize', fontSize);

subtitle([ '$K=[10^4,10^{4.5},\ldots,10^6]$ and ' ...
    '$\ell=[1.05,1.20,\ldots,1.95]$'], ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;
%xlim([1e4 1e6]);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

exportgraphics(gcf, 'energyH1evolution.pdf', 'ContentType', 'vector');

%% ========================================================================
%  3. H^2 ENERGY EVOLUTION
% ========================================================================

figure('Position', figPosition);
hold on;

for j = numel(K_values):-1:1
    for i = numel(ells):-1:1
    
        ell = ells(i);
        Kstr = K_file_strings{j};
    
        filename = sprintf('*K_%s_ell1_%.2f_ell2_%.2f_*.dat', Kstr, ell, ell);
    
        try
            data = readmatrix(dir(filename).name);
            timept = data(:,1);
            energy = data(:,4);
            loglog( timept, energy, '-', ...
                    'Color', colors(i,:), ...
                    'LineWidth', 1.5, ...
                    'MarkerSize', 6);
        catch
        end
    end
end
xlabel('Time $t$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

ylabel('$H^2$ Energy $\| \phi(t;\widetilde \varphi) \|^2_{H^2}$', ...
       'Interpreter', 'latex', ...
       'FontSize', fontSize);

title('$H^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky', ...
      'Interpreter', 'latex', ...
      'FontSize', fontSize);

subtitle([ '$K=[10^4,10^{4.5},\ldots,10^6]$ and ' ...
    '$\ell=[1.05,1.20,\ldots,1.95]$'], ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;
%xlim([1e4 1e6]);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

exportgraphics(gcf, 'energyH2evolution.pdf', 'ContentType', 'vector');
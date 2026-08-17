clear;
close all;
clc;

% MATLAB replacement for:
%   plot.gp
%   plotbranch.txt
%   plotK.txt
%   plotL.txt
%   plotTK.txt
%   plotTL.txt
%
% Run this script from the directory containing the .dat files.
%
% Outputs:
%   powerlawK.pdf
%   powerlawL.pdf
%   powerlawTK.pdf
%   powerlawTL.pdf
%   branch.pdf

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
%  1. ENERGY SCALING WITH INITIAL ENERGY K
%     gnuplot: plotK.txt
%
%     Fit:
%         log10(K*) = b + a log10(K)
%     so:
%         K* ~ K^a
% ========================================================================

figure('Position', figPosition);
hold on;

aK = zeros(size(ells));
bK = zeros(size(ells));

for i = 1:numel(ells)

    ell = ells(i);

    filename = sprintf( ...
        'powerlawK_IC_s1_ell1_%.2f_ell2_%.2f.dat', ...
        ell, ell);

    data = readmatrix(filename);

    K = data(:,1);
    Kstar = data(:,7);

    % Remove invalid entries before taking logarithms.
    valid = isfinite(K) & isfinite(Kstar) & K > 0 & Kstar > 0;
    Kfit = K(valid);
    KstarFit = Kstar(valid);

    % Equivalent to:
    % fit f(x) ... using (log10($1)):(log10($7))
    p = polyfit(log10(Kfit), log10(KstarFit), 1);

    aK(i) = p(1);
    bK(i) = p(2);

    loglog( ...
        K, Kstar, '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf( ...
            '$\\widetilde K_{%.2f} \\approx 10^{%.3f} K^{%.3f}$ \\,', ...
            ell, bK(i), aK(i)));
end

powerlawK = mean(aK);
stdevK = std(aK, 0);   % sample standard deviation, denominator N-1

xlabel('Initial Energy $K$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

ylabel('Peak Transient Energy $\widetilde{K}_{\ell}$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

title('$L^2$ Energy-Scaling Law for 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

subtitle(sprintf([ ...
    '$K=[10^4,10^{4.5},\\ldots,10^6]$ and ' ...
    '$\\ell=[1.05,1.20,\\ldots,1.95]$: ' ...
    '$\\widetilde{K}_{\\ell} \\sim K^{%.4f \\pm %.4f}$'], ...
    powerlawK, stdevK), ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;
xlim([1e4 1e6]);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

legend( ...
    'Location', 'southeast', ...
    'FontSize', fontSize-4, ...
    'Interpreter', 'latex');

exportgraphics(gcf, 'powerlawK.pdf', 'ContentType', 'vector');

fprintf('Energy scaling K exponents:\n');
disp(aK);
fprintf('Mean exponent = %.8f\n', powerlawK);
fprintf('Sample std    = %.8f\n\n', stdevK);


%% ========================================================================
%  2. ENERGY SCALING WITH DOMAIN FACTOR ell
%     gnuplot: plotL.txt
%
%     Fit:
%         log10(K*) = b + a ell
%     so:
%         K* ~ 10^(a ell + b)
% ========================================================================

figure('Position', figPosition);
hold on;

aL = zeros(size(K_values));
bL = zeros(size(K_values));

for i = 1:numel(K_values)

    filename = sprintf( ...
        'powerlawL_IC_s1_K_%s.dat', ...
        K_file_strings{i});

    data = readmatrix(filename);

    ell = data(:,2);
    Kstar = data(:,7);

    valid = isfinite(ell) & isfinite(Kstar) & ell > 0 & Kstar > 0;
    ellFit = ell(valid);
    KstarFit = Kstar(valid);

    % Equivalent to:
    % fit f(x) ... using 2:(log10($7))
    p = polyfit(ellFit, log10(KstarFit), 1);

    aL(i) = p(1);
    bL(i) = p(2);

    semilogy( ...
        ell, Kstar, '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf( ...
            '$\\widetilde K_{%s} \\approx 10^{%.3f}10^{%.3f\\ell}$', ...
            K_labels{i}, bL(i), aL(i)));
end

powerlawL = mean(aL);
stdevL = std(aL, 0);

xlabel('Domain Factor $\ell$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

ylabel('Peak Transient Energy $\widetilde{K}_{K}$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

title('Isotropic Domain-Scaling Law for 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

subtitle(sprintf([ ...
    '$K=[10^4,10^{4.5},\\ldots,10^6]$ and ' ...
    '$\\ell=[1.05,1.20,\\ldots,1.95]$: ' ...
    '$\\widetilde{K}_{K} \\sim 10^{%.4f\\ell \\pm %.4f}$'], ...
    powerlawL, stdevL), ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box off;
xlim([1 2]);

set(gca, ...
    'XScale', 'linear', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

legend( ...
    'Location', 'southoutside', ...
    'FontSize', fontSize-4, ...
    'Interpreter', 'latex',...
    'NumColumns', 3);

exportgraphics(gcf, 'powerlawL.pdf', 'ContentType', 'vector');

fprintf('Domain scaling energy exponents:\n');
disp(aL);
fprintf('Mean exponent = %.8f\n', powerlawL);
fprintf('Sample std    = %.8f\n\n', stdevL);


%% ========================================================================
%  3. TIME SCALING WITH INITIAL ENERGY K
%     gnuplot: plotTK.txt
%
%     Fit:
%         log10(T*) = b + a log10(K)
%     so:
%         T* ~ K^a
% ========================================================================

figure('Position', figPosition);
hold on;

aTK = zeros(size(ells));
bTK = zeros(size(ells));

for i = 1:numel(ells)

    ell = ells(i);

    filename = sprintf( ...
        'powerlawTK_IC_s1_ell1_%.2f_ell2_%.2f.dat', ...
        ell, ell);

    data = readmatrix(filename);

    K = data(:,1);
    Tstar = data(:,4);

    valid = isfinite(K) & isfinite(Tstar) & K > 0 & Tstar > 0;
    Kfit = K(valid);
    TstarFit = Tstar(valid);

    % Equivalent to:
    % fit f(x) ... using (log10($1)):(log10($4))
    p = polyfit(log10(Kfit), log10(TstarFit), 1);

    aTK(i) = p(1);
    bTK(i) = p(2);

    loglog( ...
        K, Tstar, '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf( ...
            '$\\widetilde{T}_{%.2f} \\approx 10^{%.3f}K^{%.3f}$', ...
            ell, bTK(i), aTK(i)));
end

powerlawTK = mean(aTK);
stdevTK = std(aTK, 0);

xlabel('Initial Energy $K$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

ylabel('Time Window of Peak Transient Energy $\widetilde{T}_{\ell}$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

title('$L^2$ Energy Time-Scaling Law for 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

subtitle(sprintf([ ...
    '$K=[10^4,10^{4.5},\\ldots,10^6]$ and ' ...
    '$\\ell=[1.05,1.20,\\ldots,1.95]$: ' ...
    '$\\widetilde{T}_{\\ell} \\sim K^{%.4f \\pm %.4f}$'], ...
    powerlawTK, stdevTK), ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;
xlim([1e4 1e6]);

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

legend( ...
    'Location', 'southwest', ...
    'FontSize', fontSize-4, ...
    'Interpreter', 'latex');

exportgraphics(gcf, 'powerlawTK.pdf', 'ContentType', 'vector');

fprintf('Time-vs-K exponents:\n');
disp(aTK);
fprintf('Mean exponent = %.8f\n', powerlawTK);
fprintf('Sample std    = %.8f\n\n', stdevTK);


%% ========================================================================
%  4. TIME SCALING WITH DOMAIN FACTOR ell
%     gnuplot: plotTL.txt
%
%     Fit:
%         log10(T*) = b + a ell
%     so:
%         T* ~ 10^(a ell + b)
% ========================================================================

figure('Position', figPosition);
hold on;

aTL = zeros(size(K_values));
bTL = zeros(size(K_values));

for i = 1:numel(K_values)

    filename = sprintf( ...
        'powerlawTL_IC_s1_K_%s.dat', ...
        K_file_strings{i});

    data = readmatrix(filename);

    ell = data(:,2);
    Tstar = data(:,4);

    valid = isfinite(ell) & isfinite(Tstar) & ell > 0 & Tstar > 0;
    ellFit = ell(valid);
    TstarFit = Tstar(valid);

    % Equivalent to:
    % fit f(x) ... using 2:(log10($4))
    p = polyfit(ellFit, log10(TstarFit), 1);

    aTL(i) = p(1);
    bTL(i) = p(2);

    semilogy( ...
        ell, Tstar, '-o', ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf( ...
            '$\\widetilde{T}_{%s} \\approx 10^{%.3f}10^{%.3f\\ell}$', ...
            K_labels{i}, bTL(i), aTL(i)));

end

powerlawTL = mean(aTL);
stdevTL = std(aTL, 0);

xlabel('Domain Factor $\ell$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

ylabel('Time Window of Peak Transient Energy $\widetilde{T}_{K}$', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

title('Isotropic Domain Time-Scaling Law for 2D Kuramoto-Sivashinsky', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

subtitle(sprintf([ ...
    '$K=[10^4,10^{4.5},\\ldots,10^6]$ and ' ...
    '$\\ell=[1.05,1.20,\\ldots,1.95]$: ' ...
    '$\\widetilde{T}_{K} \\sim 10^{%.4f\\ell \\pm %.4f}$'], ...
    powerlawTL, stdevTL), ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box off;
xlim([1 2]);

set(gca, ...
    'XScale', 'linear', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

legend( ...
    'Location', 'southeast', ...
    'FontSize', fontSize-4, ...
    'Interpreter', 'latex', ...
    'NumColumns', 3);

exportgraphics(gcf, 'powerlawTL.pdf', 'ContentType', 'vector');

fprintf('Time-vs-domain exponents:\n');
disp(aTL);
fprintf('Mean exponent = %.8f\n', powerlawTL);
fprintf('Sample std    = %.8f\n\n', stdevTL);


%% ========================================================================
%  5. SOLUTION BRANCHES
%     gnuplot: plotbranch.txt
%
%     Column 4 = time window
%     Column 7 = maximum L2 energy amplification
%
%     The original gnuplot file plots all 5 K values x all 7 ell values.
%     Curves with the same ell use the same color.
% ========================================================================

figure('Position', figPosition);
hold on;

legendHandles = gobjects(numel(ells),1);

for iK = 1:numel(K_values)

    for iEll = 1:numel(ells)

        ell = ells(iEll);

        filename = sprintf( ...
            'branch_IC_s1_K_%s_ell1_%.2f_ell2_%.2f.dat', ...
            K_file_strings{iK}, ell, ell);

        data = readmatrix(filename);

        T = data(:,4);
        Kstar = data(:,7);

        h = loglog( ...
            T, Kstar, '-o', ...
            'Color', colors(iEll,:), ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5);

        % Only one legend entry per ell, matching the conceptual intent
        % of the gnuplot plot where color denotes ell.
        if iK == 1
            h.DisplayName = sprintf('%.2f', ell);
            legendHandles(iEll) = h;
        else
            h.HandleVisibility = 'off';
        end
    end
end

xlabel('Time Window', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

ylabel('Maximum $L^2$ Energy Amplification', ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

title([ ...
    'Solution Branches for $K=[10^4,10^{4.5},\ldots,10^6]$ and ' ...
    '$\ell=[1.05,1.20,\ldots,1.95]$'], ...
    'Interpreter', 'latex', ...
    'FontSize', fontSize);

grid on;
box on;

set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', fontSize);

legend( ...
    legendHandles, ...
    'Location', 'northeast', ...
    'NumColumns', 5, ...
    'FontSize', fontSize, ...
    'Interpreter', 'latex');

exportgraphics(gcf, 'branch.pdf', 'ContentType', 'vector');


%% ========================================================================
%  Summary
% ========================================================================

fprintf('\n============================================================\n');
fprintf('SUMMARY OF FITTED SCALINGS\n');
fprintf('============================================================\n');

fprintf('K* vs K:    exponent = %.8f +/- %.8f\n', ...
    powerlawK, stdevK);

fprintf('K* vs ell:  log10 slope = %.8f +/- %.8f\n', ...
    powerlawL, stdevL);

fprintf('T* vs K:    exponent = %.8f +/- %.8f\n', ...
    powerlawTK, stdevTK);

fprintf('T* vs ell:  log10 slope = %.8f +/- %.8f\n', ...
    powerlawTL, stdevTL);

fprintf('============================================================\n');

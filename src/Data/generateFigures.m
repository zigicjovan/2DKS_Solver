close all

testcase = '_IC_s1_N1_64_N2_64_dt_1.0e-04_K_3.2e+00_ell1_1.02_ell2_1.02_T_3.16e-02_opt_1_tol_1e-06_cont_1_optT_3.16e-02';
forwardDir = 'ForwardSolution';
forwardFile = 'fwd_';
energyDir = 'EnergyEvolution';
energyFile = 'energy';
spectrumDir = 'FourierSpectrumEvolution';
spectrumFile = 'spectrum';
initialconditionDir = 'InitialData';
initialconditionFile = 'fwdIC';
terminalconditionDir = 'TerminalData';
terminalconditionFile = 'fwdTC';

root = pwd;
addpath(genpath(fullfile(root,testcase)));

numberOfStates = 100;

% Match parameter names followed by numerical values
tokens = regexp(testcase, ...
    '(?<name>N1|N2|dt|K|ell1|ell2|T|opt|tol|cont|optT)_(?<value>[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', ...
    'names');

% Store the values in a structure
params = struct();

for i = 1:numel(tokens)
    params.(tokens(i).name) = str2double(tokens(i).value);
end

gridSize = params.N1 * params.N2;

%% Find and numerically sort files by the final number
forwardFiles = dir(fullfile([testcase '/' forwardDir], [forwardFile, '*.dat']));
fileTimes = zeros(numel(forwardFiles), 1);
for k = 1:numel(forwardFiles)
    token = regexp(forwardFiles(k).name, '_([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\.dat$', 'tokens', 'once');
    if isempty(token)
        error('Could not extract final number from: %s', forwardFiles(k).name);
    end
    fileTimes(k) = str2double(token{1});
end
[fileTimes, order] = sort(fileTimes);
forwardFiles = forwardFiles(order);

energyFiles = dir(fullfile([testcase '/' energyDir], [energyFile, '*.dat']));
spectrumFiles = dir(fullfile([testcase '/' spectrumDir], [spectrumFile, '*.dat']));
initialconditionFiles = dir(fullfile([testcase '/' initialconditionDir], [initialconditionFile, '*.dat']));
terminalconditionFiles = dir(fullfile([testcase '/' terminalconditionDir], [terminalconditionFile, '*.dat']));

%% plot settings and handles
wordsize = 16;
number_of_plots = 7;
maxplotcols = 4;
number_of_plot_rows = ceil(number_of_plots/maxplotcols);
figwidth = 400*maxplotcols;
figheight = 450*number_of_plot_rows;

fig = figure('Visible', 'off');
set(fig, 'Position', [100 100 figwidth figheight], 'Color', 'white', 'Resize', 'off');

set(groot, ...
'DefaultAxesLooseInset',            [0 0 0 0], ...
'DefaultLegendFontSize',            ceil(3*wordsize/4), ...
'DefaultTextInterpreter',           'latex', ...
'DefaultAxesTickLabelInterpreter',  'latex', ...
'DefaultLegendInterpreter',         'latex');

h_surf1 = []; h_surf2 = []; a_xline = []; 
h_spec1 = []; h_proj = []; h_spec2 = [];
proj_strip = []; proj_coeffpts = [];
didInitPlots = false;

%% set field coordinates
xField = linspace(0, params.ell1, params.N1);
yField = linspace(0, params.ell2, params.N2);

% 4-period physical space domain
x2_pts = linspace( 0 , 2*params.ell1 , 2*params.N1 ); 
y2_pts = linspace( 0 , 2*params.ell2 , 2*params.N2 ); 
[ x2x , y2x ] = meshgrid(x2_pts,y2_pts); % 2-dimensional grid

%% set energy coordinates
filename = fullfile(energyFiles(k).folder, energyFiles(k).name);
fid = fopen(filename, 'r');
energyData = fscanf(fid, '%f', [4, Inf]).';
fclose(fid);
numPoints = size(energyData, 1);
xEnergy = linspace(0, energyData(end,1), numPoints);
ymin_energy = 0.5 * min(min(energyData(:,2:4)));
ymax_energy = 1.5 * max(max(energyData(:,2:4)));

energyIndices = round(linspace(1, size(energyData, 1), numberOfStates));
energyDataSampled = energyData(energyIndices, :);

%% set spectrum coordinates
k1 = [0:params.N1/2-1, -params.N1/2:-1] / params.ell1;
k2 = [0:params.N2/2-1, -params.N2/2:-1] / params.ell2;
maxRadius = hypot(max(abs(k1)), max(abs(k2)));
numRadialBins = round(maxRadius) + 1;

filename = fullfile(spectrumFiles(k).folder, spectrumFiles(k).name);
fid = fopen(filename, 'r');
spectrumData = fscanf(fid, '%f', [numberOfStates + 1, Inf]).';
spectrumData = spectrumData(2:end,:);
fclose(fid);
xmax_spec = size(spectrumData,1);
ymax_spec = 1.5 * max(max(spectrumData));

%% create spectrum IC figure
currentFig = figure;
set(currentFig, 'Position', [100 100 900 750]);
set(currentFig, 'Color', 'white');

ax = axes(currentFig);
set(ax, 'FontSize', wordsize);
set(ax, 'Color', 'white');

semilogy(ax(k), 1:xmax_spec, spectrumData(:,2), "--");
xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$' ); 
legend(ax(k),'$E(k)$','Box','off','FontSize',wordsize);
title(ax(k),"Energy spectrum" );
xlim(ax(k), [1 xmax_spec]);
ylim(ax(k), [1e-20 ymax_spec]);
axis(ax(k),'square');

filename = fullfile(spectrumDir, ['spectrumIC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

%% create spectrum TC figure
currentFig = figure;
set(currentFig, 'Position', [100 100 900 750]);
set(currentFig, 'Color', 'white');

ax = axes(currentFig);
set(ax, 'FontSize', wordsize);
set(ax, 'Color', 'white');

semilogy(ax(k), 1:xmax_spec, spectrumData(:,end), "--");
xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$' ); 
%legend(ax(k),'$E(k)$','$Ce^{-2\delta k}$','Box','off','FontSize',wordsize);
legend(ax(k),'$E(k)$','Box','off','FontSize',wordsize);
title(ax(k),"Energy spectrum" );
xlim(ax(k), [1 xmax_spec]);
ylim(ax(k), [1e-20 ymax_spec]);
axis(ax(k),'square');

filename = fullfile(spectrumDir, ['spectrumTC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

%% AS strip width fit
astripwidth = NaN(numberOfStates,2);
for i = 1:numberOfStates    
    % Width of the analyticity strip [C , delta]
    if i > 1
        [astripwidth(i,1),astripwidth(i,2)] = expfit_delta(spectrumData(2:end,1), spectrumData(2:end,i+1), astripwidth(i-1,2));
    else
        [astripwidth(i,1),astripwidth(i,2)] = expfit_delta(spectrumData(2:end,1), spectrumData(2:end,i+1), 1e15);
    end
    asstrip_fit = spectrumData(:,2:end);
end

%% process analyticity strip width
for i = 2:size(spectrumData,2)
    asstrip_fit(:,i-1) = astripwidth(i-1,1)*exp(-2*astripwidth(i-1,2)*spectrumData(:,1));
end

%% create field IC figure
currentFig = figure;
set(currentFig, 'Position', [100 100 900 750]);
set(currentFig, 'Color', 'white');

ax = axes(currentFig);
set(ax, 'FontSize', wordsize);
set(ax, 'Color', 'white');

filename = fullfile(initialconditionFiles(k).folder, initialconditionFiles(k).name);
fid = fopen(filename, 'rb');
raw = fread(fid, Inf, 'double');
fclose(fid);

i = 1;
initialField = raw(1:2:end) + 1i * raw(2:2:end);
states = reshape(initialField, gridSize, []);
states(:, all(states == 0, 1)) = [];
initialField = states(:);
u_hat = reshape(initialField(1:gridSize), params.N1, params.N2);
u = real(ifft2(u_hat));

%% appropriate phase-shift to center any local minima
minctr = 1;
[~,minidx] = mink(u(:),20);
while minctr < size(minidx,1)
    minidx = [minidx(1:minctr) ; minidx(abs(minidx(minctr) - minidx(minctr:end)) > .05*size(u(:),1))];
    minctr = minctr + 1;
end
[row, col] = ind2sub(size(u), minidx);
sz = size(u);

% periodic mean for rows
theta_r = 2*pi*(row-1)/sz(1);
rowmin = mod(atan2(mean(sin(theta_r)), mean(cos(theta_r))) * sz(1)/(2*pi), sz(1)) + 1;

% periodic mean for columns
theta_c = 2*pi*(col-1)/sz(2);
colmin = mod(atan2(mean(sin(theta_c)), mean(cos(theta_c))) * sz(2)/(2*pi), sz(2)) + 1;

rowmin = round(rowmin);
colmin = round(colmin);
[Ny, Nx] = size(u);
rowc = floor((Ny+1)/2) + 1;
colc = floor((Nx+1)/2) + 1;
drow = rowc - rowmin;
dcol = colc - colmin;
u = circshift(u, [drow, dcol]);
%%

surfc(ax, xField, yField, u);
view(ax, 2);
shading(ax, 'interp');
colormap(ax, redblue);
colorbar(ax);
xlabel(ax, '$\frac{x_1}{2\pi}$', 'Interpreter', 'latex');
ylabel(ax, '$\frac{x_2}{2\pi}$', 'Interpreter', 'latex');
axis(ax, 'equal', 'tight');

currentT = fileTimes(k)*(i-1)/(numberOfStates-1);
titleText = sprintf( ...
                ['${N_1 = %d, N_2 = %d, \\Delta t = %.1e, K = %.0e,' ...
                 ' \\ell_1 = %.2f, \\ell_2 = %.2f, T = %.6f}$'], ...
                params.N1, params.N2, params.dt, params.K, ...
                params.ell1, params.ell2, currentT);
title(ax, titleText, 'Interpreter', 'latex');

drawnow;

filename = fullfile(initialconditionDir, ['contourIC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

view(ax, 3);
axis(ax,"square")
colorbar(ax, 'off');
filename = fullfile(initialconditionDir, ['surfaceIC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

%% create field TC figure
currentFig = figure;
set(currentFig, 'Position', [100 100 900 750]);
set(currentFig, 'Color', 'white');

ax = axes(currentFig);
set(ax, 'FontSize', wordsize);
set(ax, 'Color', 'white');

filename = fullfile(terminalconditionFiles(k).folder, terminalconditionFiles(k).name);
fid = fopen(filename, 'rb');
raw = fread(fid, Inf, 'double');
fclose(fid);

i = numberOfStates;
terminalField = raw(1:2:end) + 1i * raw(2:2:end);
states = reshape(terminalField, gridSize, []);
states(:, all(states == 0, 1)) = [];
terminalField = states(:);

indexStart = (i - 1) * gridSize + 1;
indexEnd   = i * gridSize;

u_hat = reshape(terminalField(1:gridSize), params.N1, params.N2);
u = real(ifft2(u_hat));
u = circshift(u, [drow, dcol]);

surfc(ax, u);
view(ax, 2);
shading(ax, 'interp');
colormap(ax, redblue);
colorbar(ax);
axis(ax, 'equal', 'tight');

currentT = fileTimes(k)*(i-1)/(numberOfStates-1);
titleText = sprintf( ...
                ['${N_1 = %d, N_2 = %d, \\Delta t = %.1e, K = %.0e,' ...
                 ' \\ell_1 = %.2f, \\ell_2 = %.2f, T = %.6f}$'], ...
                params.N1, params.N2, params.dt, params.K, ...
                params.ell1, params.ell2, currentT);
title(ax, titleText, 'Interpreter', 'latex');

drawnow;

filename = fullfile(terminalconditionDir, ['contourTC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

view(ax, 3);
axis(ax,"square")
colorbar(ax, 'off');
filename = fullfile(terminalconditionDir, ['surfaceTC' testcase]);
saveas(currentFig,[filename '.fig'])
exportgraphics(currentFig,[filename '.pdf'])

%% Create movie
gifFile = fullfile(forwardDir, ['movie' testcase '.gif']);
frameDelay = 1/20;  % equivalent to 20 frames per second
isFirstFrame = true;

%% Precompute projection coefficients for all states

% Copy only the radial-bin radii, not the spectrum values
numberOfStates = size(spectrumData,2) - 1;
projcoeffradialevolution = [spectrumData(:,1), zeros(size(spectrumData,1),numberOfStates)];

radialModeLabels = strings(size(projcoeffradialevolution,1),1);
projcoeffmodeevolution = [];

stateIndex = 0;

for fi = 1:numel(forwardFiles)
    filename = fullfile(forwardFiles(fi).folder,forwardFiles(fi).name);

    fid = fopen(filename,'rb');
    raw = fread(fid,Inf,'double');
    fclose(fid);

    fwdField = raw(1:2:end) + 1i*raw(2:2:end);
    fileStates = reshape(fwdField,gridSize,[]);

    % Remove unused preallocated zero states
    fileStates(:,all(fileStates == 0,1)) = [];

    for localIndex = 1:numberOfStates
        stateIndex = stateIndex + 1;

        u_hat = reshape(fileStates(:,localIndex),params.N1,params.N2);

        u = real(ifft2(u_hat));
        u = circshift(u,[drow,dcol]);

        [~,projcoeffs,unstablemodes] = eigenfunction_validation(u(:),params.ell1,params.N1,'full');

        projectionWeights = abs(projcoeffs).^2;
        weightSum = sum(projectionWeights);

        if weightSum > 0
            projectionWeights = projectionWeights/weightSum;
        else
            projectionWeights(:) = 0;
        end

        if stateIndex == 1
            numberOfModes = size(unstablemodes,1);
            projcoeffmodeevolution = zeros(numberOfModes,numberOfStates + 2);
            projcoeffmodeevolution(:,1:2) = unstablemodes(:,1:2);

            % Exact physical radius of every unstable mode
            modeRadii = hypot(unstablemodes(:,1)/params.ell1, unstablemodes(:,2)/params.ell2);
        
            % Group sign changes and symmetric modes having equal radii
            [radialRadii,representativeIdx,modeToRadialGroup] = uniquetol(modeRadii,1e-10,'DataScale',1);
        
            numberOfRadialGroups = numel(radialRadii);
            projcoeffradialevolution = [radialRadii, zeros(numberOfRadialGroups,numberOfStates)];
            representativeModes = unstablemodes(representativeIdx,1:2);
        
            % Make representative labels consistent:
            % axial modes as (0,n), mixed modes with larger index first.
            representativeModes = abs(representativeModes);
        
            for j = 1:size(representativeModes,1)
                if all(representativeModes(j,:) > 0) && representativeModes(j,1) < representativeModes(j,2)
                    representativeModes(j,:) = fliplr(representativeModes(j,:));
                elseif representativeModes(j,2) == 0
                    representativeModes(j,:) = fliplr(representativeModes(j,:));
                end
            end
        
            modelabels = "(" + string(representativeModes(:,1)) + "," + string(representativeModes(:,2)) + ")"; 
            modecats = categorical(modelabels,modelabels,'Ordinal',true);
        elseif ~isequal(unstablemodes(:,1:2), projcoeffmodeevolution(:,1:2))
            error('Unstable-mode ordering changed at state %d.', stateIndex);
        end

        projcoeffmodeevolution(:,stateIndex+2) = projectionWeights;

        % Sum equivalent modes into their exact radial groups
        radialWeights = accumarray(modeToRadialGroup, projectionWeights, [numberOfRadialGroups,1], @sum,0);        
        projcoeffradialevolution(:,stateIndex+1) = radialWeights;
    end
end

projcoeffradialcut = projcoeffradialevolution;

% Choose fixed limits appropriate for your solution.
close all
colorLimits = [-100, 100];

currentFig = figure;
set(currentFig, 'Color', 'white');

for k = 1:number_of_plots
    ax(k) = subplot(number_of_plot_rows,maxplotcols,k);
end

for k = 1:number_of_plots
    set(ax(k),'Color','white');
    axis(ax(k), 'square');
    pos = get(ax(k),'Position');
    pos(2) = pos(2) - 0.05;
    set(ax(k),'Position',pos);
end

title1 = 'Forward-time 2DKS solution';
title2 = '';
sg = sgtitle({title1, title2}, 'Interpreter','latex','FontSize',ceil(1.25*wordsize));

stateIndex = 0;

for fi = 1:numel(forwardFiles)

    filename = fullfile(forwardFiles(fi).folder, forwardFiles(fi).name);
    fid = fopen(filename, 'rb');
    raw = fread(fid, Inf, 'double');
    fclose(fid);

    % Convert alternating real and imaginary values to complex values
    fwdField = raw(1:2:end) + 1i * raw(2:2:end);
    states = reshape(fwdField, gridSize, []);
    states(:, all(states == 0, 1)) = [];
    fwdField = states(:);

    % Number of complete states contained in this file
    for i = 1:numberOfStates
        stateIndex = stateIndex + 1;
        currentT = energyDataSampled(i,1);

        %% prepare physical field
        indexStart = (i - 1) * gridSize + 1;
        indexEnd   = i * gridSize;

        u_hat = reshape(fwdField(indexStart:indexEnd), params.N1, params.N2);
        u = real(ifft2(u_hat));
        u = circshift(u, [drow, dcol]);
        u_2x = [u , u ; u , u];
        
        %% produce subplots
        if stateIndex == 1
            % ------------ AXIS 1: physical field ---------------
            k = 1;
            h_surf1 = surfc(ax(k), xField, yField, u);
            xlabel(ax(k),'$\frac{x_1}{2\pi}$');
            ylabel(ax(k),'$\frac{x_2}{2\pi}$');
            %title(ax(k),"Periodic solution field");
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            view(ax(k),3);
            pbaspect(ax(k), [ params.ell1 params.ell2 max(params.ell1, params.ell2) ])
        
            % ------------ AXIS 2: tiled domain -----------------
            k = 2;
            h_surf2 = surfc(ax(k), x2x, y2x, u_2x);
            xlabel(ax(k),'$\frac{x_1}{2\pi}$' );
            ylabel(ax(k),'$\frac{x_2}{2\pi}$' );
            title(ax(k),"Tiled solution field" );
            shading(ax(k),'interp');
            colormap(ax(k), redblue);
            xline(ax(k), params.ell1, '--');
            yline(ax(k), params.ell2, '--');
            view(ax(k),2);
            axis(ax(k), [0 2*params.ell1 0 2*params.ell2])
            axis(ax(k), 'equal')

            % ------------ AXIS 3: energy evolution -------------
            k = 3;
            semilogy(ax(k), energyData(:,1), energyData(:,4), 'g'); 
            hold(ax(k),'on');
            semilogy(ax(k), energyData(:,1), energyData(:,3), 'r');
            semilogy(ax(k), energyData(:,1), energyData(:,2), 'b');
            H2_xline = plot(ax(k), energyData(i,1), energyData(i,4), 'ko');
            H1_xline = plot(ax(k), energyData(i,1), energyData(i,3), 'ko');
            L2_xline = plot(ax(k), energyData(i,1), energyData(i,2), 'ko');
            hold(ax(k),'off');
            
            xlabel(ax(k),'Time $t$' );
            ylabel(ax(k),'$\| \phi(t;\varphi) \|^2_{S}$' );
            xlim(ax(k), [0 energyData(end,1)]);
            ylim(ax(k), [ymin_energy ymax_energy]);
            title(ax(k), "Energy");
            legend(ax(k), '$S=H^2$','$S=H^1$','$S=L^2$','Location','southeast','Box','off');
            axis(ax(k),'square');

            % ------------ AXIS 4: radial spectrum --------------
            k = 4;
            h_spec1 = semilogy(ax(k), 1:xmax_spec, spectrumData(:,i+1), "--");
            hold(ax(k),'on');
            h_spec2 = semilogy(ax(k), 1:xmax_spec, asstrip_fit(:,i), "r--");
            hold(ax(k),'off');
            xlabel(ax(k),'$k \approx \sqrt{k_1^2+k^2_2}$' ); 
            legend(ax(k),'$E(k)$','$Ce^{-2\delta k}$','Box','off','FontSize',wordsize);
            title(ax(k),"Spectrum" );
            xlim(ax(k), [1 xmax_spec]);
            ylim(ax(k), [1e-20 ymax_spec]);
            axis(ax(k),'square');

            % ------------ AXIS 5: analyticity strip -----------
            k = 5;
            plot(ax(k), energyDataSampled(:,1), astripwidth(:,2), 'b'); 
            hold(ax(k),'on');
            a_xline = plot(ax(k), energyDataSampled(stateIndex,1), astripwidth(stateIndex,2), 'ko');
            yline(ax(k), max(2*pi*params.ell1/params.N1,2*pi*params.ell2/params.N2), '--');
            ylim(ax(k), [0 1.5*max(astripwidth(:,2))]);
            xlim(ax(k), [0 energyDataSampled(end,1)]);
            hold(ax(k),'off');
            
            xlabel(ax(k),'Time $t$' );
            ylabel(ax(k),'$\delta(t)$' );
            title(ax(k), "Analyticity Strip Width");
            axis(ax(k),'square');

            % ------------ AXIS 6: projection coefficients by radial group ------------
            k = 6;
            cla(ax(k));
            
            numberOfRadialGroups = size(projcoeffradialcut,1);
            xmodelabels = 1:numberOfRadialGroups;
            
            h_proj = bar(ax(k), xmodelabels, projcoeffradialcut(:,stateIndex+1), 'BarWidth',1.0);
            
            ax(k).XTick = xmodelabels;
            ax(k).XTickLabel = cellstr(modecats);
            ax(k).XLim = [0.5,numberOfRadialGroups + 0.5];        
            xlabel(ax(k),'$(k_1,k_2)$', 'Interpreter','latex');
            ylabel(ax(k),'$P(a_k)$', 'Interpreter','latex');            
            title(ax(k),'Projection coefficient weights');
            
            maxRadialWeight = max(projcoeffradialcut(:,2:end),[],'all');
            
            if maxRadialWeight > 0
                ylim(ax(k),[0,1.1*maxRadialWeight]);
            else
                ylim(ax(k),[0,1]);
            end
            
            axis(ax(k),'square');

            % ------------ AXIS 7: radial projection coefficients v time ------------
            k = 7;
            proj_strip = plot(ax(k), energyDataSampled(:,1), projcoeffradialcut(:,2:end)');            
            maxRadialWeight = max(projcoeffradialcut(:,2:end),[],'all');
            
            if maxRadialWeight > 0
                modalfigylim = [0,1.05*maxRadialWeight];
            else
                modalfigylim = [0,1];
            end
            
            hold(ax(k),'on');
            numberOfRadialGroups = size(projcoeffradialcut,1);
            projcoeff_xpts = currentT*ones(numberOfRadialGroups,1);
            projcoeff_ypts = projcoeffradialcut(:,stateIndex+1);
            proj_coeffpts = plot(ax(k), projcoeff_xpts, projcoeff_ypts, 'ko', 'LineStyle','none');
            
            % One label for each radial group
            labelsToPlot = modelabels;
            labelY = projcoeffradialcut(:,end);
            
            [labelY,order] = sort(labelY);
            labelsToPlot = labelsToPlot(order);
            
            figfrac = 0.07*diff(modalfigylim);
            
            for imode = 2:numel(labelY)
                if labelY(imode) - labelY(imode-1) < figfrac
                    labelY(imode) = labelY(imode-1) + figfrac;
                end
            end
            
            modalfigylim(2) = max(modalfigylim(2), max(labelY) + 0.5*figfrac);
            
            % Keep the axes box ending at final T
            xlim(ax(k),[0,energyDataSampled(end,1)]);
            
            for imode = 1:numel(labelsToPlot)
                text(ax(k),energyDataSampled(end,1),labelY(imode),labelsToPlot(imode), ...
                    'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                    'Interpreter','latex', 'FontSize',ceil(0.75*wordsize), 'Clipping','off');
            end
            
            xlim(ax(k),[0,energyDataSampled(end,1)]);
            ylim(ax(k),modalfigylim);
            
            hold(ax(k),'off');
            
            xlabel(ax(k),'Time $t$');
            ylabel(ax(k),'$P(a_k)$');
            title(ax(k),'Radial modal energy');
            axis(ax(k),'square');

            % ---- enforce fonts/interpreters AFTER all plotting calls ----
            for kk = 1:number_of_plots
            set(ax(kk), 'FontSize', ceil(0.75*wordsize), 'LabelFontSizeMultiplier', 1.0, ...
                'TitleFontSizeMultiplier', 1.0, 'TickLabelInterpreter', 'latex');  % axes tick labels
            
            % Ensure existing labels/titles use LaTeX + correct size
            ax(kk).XLabel.Interpreter = 'latex';
            ax(kk).YLabel.Interpreter = 'latex';
            ax(kk).Title.Interpreter  = 'latex';
            
            ax(kk).XLabel.FontSize = wordsize;
            ax(kk).YLabel.FontSize = wordsize;
            ax(kk).Title.FontSize  = wordsize;
            end
            
            % Legends: force font size/interpreter (since legend does not inherit axes)
            lg = findall(fig, 'Type', 'Legend');
            set(lg, 'Interpreter','latex', 'FontSize', ceil(0.75*wordsize));

        else
            % Axis 1 & 2 surfaces:
            set(h_surf1, 'ZData', u);
            set(h_surf2, 'ZData', u_2x);
            
            % Axis 3: just move the xline (curves are static arrays)
            set(H2_xline, 'XData', energyDataSampled(stateIndex,1), 'YData', energyDataSampled(stateIndex,4));
            set(H1_xline, 'XData', energyDataSampled(stateIndex,1), 'YData', energyDataSampled(stateIndex,3));
            set(L2_xline, 'XData', energyDataSampled(stateIndex,1), 'YData', energyDataSampled(stateIndex,2));
            
            % Axis 4: update spectrum
            set(h_spec1, 'YData', spectrumData(:,stateIndex+1));
            set(h_spec2, 'YData', asstrip_fit(:,stateIndex));
                
            % Axis 5: update analyticity strip
            set(a_xline, 'XData', energyDataSampled(stateIndex,1), 'YData', astripwidth(stateIndex,2));
            
            % Axis 6: update projection coefficient spectrum
            set(h_proj, 'YData',projcoeffradialcut(:,stateIndex+1));
            
            % Axis X: update proj coeff evolution
            k = 7;
            set(proj_coeffpts, 'XData',currentT*ones(numberOfRadialGroups,1), 'YData',projcoeffradialcut(:,stateIndex+1));
        end    

        titleText = sprintf( ...
                    ['${N_1 = %d, N_2 = %d, \\Delta t = %.1e, K = %.0e,' ...
                     ' \\ell_1 = %.2f, \\ell_2 = %.2f, T = %.6f}$'], ...
                    params.N1, params.N2, params.dt, params.K, ...
                    params.ell1, params.ell2, currentT);
        sg.String = {title1, titleText};

        drawnow;

        frame = getframe(currentFig);
        rgbImage = frame2im(frame);
        [indexedImage, colorMap] = rgb2ind(rgbImage, 256);
        
        if isFirstFrame
            imwrite(indexedImage, colorMap, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', frameDelay);
            isFirstFrame = false;
        else
            imwrite(indexedImage, colorMap, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
        end
    end
end
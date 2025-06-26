%% Figure 2: RT-share Manipulation and Violation Effects
% This script visualises how varying the first-stage RT-share in the relay model
% affects predicted CDFs, RSE, and Miller bound violations.

clear; clc; close all;

% Add custom function directory to path
addpath(fullfile(pwd, 'Functions'));

%% Load Fitted Parameters from Previous Analysis
% Parameters stored as a cell array: {A_dist, V_dist, weightStruct}
paramsFile = fullfile(pwd, 'FittedParameters', 'Miller82_Parameters.mat');
load(paramsFile, 'params');

A_dist = params{1};  % Auditory inverse Gaussian distribution
V_dist = params{2};  % Visual inverse Gaussian distribution
W      = params{3};  % Optimal relay weights

% Extract parameters for convenience
aMU     = A_dist.mu;
vMU     = V_dist.mu;
aLAMBDA = A_dist.lambda;
vLAMBDA = V_dist.lambda;
optW    = W.optimal_w;

%% Set up Grid for Weight Variation and RT Range
numWeights = 100;                  % Number of weight steps
weights    = linspace(0, 0.5, numWeights);
xx         = linspace(100, 700, 100);  % RT values for CDF evaluation

% Preallocate for speed
modelCDFs = nan(length(xx), numWeights);
gains     = nan(1, numWeights);
violationsTrue = nan(1, numWeights);

%% Load Empirical Quantiles for Plotting
loadedData = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller82', 'miller.xlsx'));
empData = loadedData(:,1:3);
probLevels = loadedData(:,4);

%% Compute Relay Model CDFs Across Weight Range
% Parallel loop speeds up independent CDF evaluations
parfor idx = 1:numWeights
    w = weights(idx);
    modelCDFs(:,idx) = getRelayCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, w, 1-w, w, 1-w);
end

%% Plotting: CDF Evolution with Weight
figure;
t = tiledlayout(1, 4, 'TileSpacing','compact', 'Padding','compact');
ax1 = nexttile(t, [1 2]); hold(ax1,'on');

linewidth = 2.5;
fontsize = 14;

% Unisensory CDFs
uniA = getUniCDF(xx, aMU, aLAMBDA);
uniV = getUniCDF(xx, vMU, vLAMBDA);
plot(xx, uniA, 'Color','#3AAA35', 'LineWidth',linewidth);
plot(xx, uniV, 'Color','#009FE3', 'LineWidth',linewidth);

% Relay model curves colored by weight
grey = [228 210 231]/255; black = [137 41 133]/255;
colors = [linspace(grey(1),black(1),numWeights)',...
          linspace(grey(2),black(2),numWeights)',...
          linspace(grey(3),black(3),numWeights)'];
for idx = 1:numWeights
    plot(ax1, xx, modelCDFs(:,idx), 'Color', colors(idx,:), 'LineWidth',1.5);
end

% Reference bounds: Miller and Raab
plot(ax1, xx, getMillerCDF(xx, aMU,vMU,aLAMBDA,vLAMBDA), 'r-', 'LineWidth',linewidth);
plot(ax1, xx, getRaabCDF(xx, aMU,vMU,aLAMBDA,vLAMBDA), 'r--', 'LineWidth',linewidth);

% Optimal-weight relay curve
plot(ax1, xx, getRelayCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA,optW,1-optW,optW,1-optW), 'k-', 'LineWidth',linewidth);

% Overlay empirical quantiles
scatter(ax1, empData(:,1), probLevels, 50, 'MarkerEdgeColor','#3AAA35', 'MarkerFaceColor','w', 'LineWidth', linewidth);
scatter(ax1, empData(:,2), probLevels, 50, 'MarkerEdgeColor','#009FE3', 'MarkerFaceColor','w', 'LineWidth', linewidth);
scatter(ax1, empData(:,3), probLevels, 50, 'MarkerEdgeColor','k',        'MarkerFaceColor','w', 'LineWidth', linewidth);

xlabel('Response Time (ms)'); ylabel('Cumulative Probability');
axis([100 700 0 1]); box off;
ax1.TickDir = 'out'; ax1.FontSize = fontsize;

% Set axis ticks and other properties
yticks(linspace(0, 1, 3))
xticks(linspace(100, 700, 3))

ax1.YColor = "k";
ax1.XColor = "k";
ax1.LineWidth = linewidth;

%% Compute RSE and Violations Across Weights
parfor idx = 1:numWeights
    w = weights(idx);
    % Race and bound CDFs
    griceCDF   = getGriceCDF(uniA, uniV);
    raabCDF    = getRaabCDF(xx, aMU,vMU,aLAMBDA,vLAMBDA);
    theCDF     = modelCDFs(:,idx);
    millerCDF  = getMillerCDF(xx, aMU,vMU,aLAMBDA,vLAMBDA);

    % Redundancy gain (RSE) and violation metric
    gains(idx)         = getRSE_fromCDF(xx, theCDF, griceCDF);
    violationsTrue(idx)= getViolation_fromCDF(xx, theCDF, millerCDF);
end

%% Plot RSE and Violation as Function of Weight
% RSE subplot
ax2 = nexttile(t); hold(ax2,'on');
plot(weights*100, gains, 'k-', 'LineWidth',linewidth);
ylabel('RSE (ms)');
xlabel('RT Share (%)');
plot([0 weights(end)*100], [getRSE_fromCDF(xx,getRaabCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), getGriceCDF(uniA,uniV))]*[1 1], 'r--','LineWidth',linewidth);
axis([0 50 0 80]); box off; ax2.TickDir='out'; ax2.FontSize=fontsize;
yticks(linspace(0, 80, 5))
xticks(linspace(0, 50, 6))
ax2.YColor = "k";
ax2.XColor = "k";
ax2.LineWidth = linewidth;


% Violation subplot
ax3 = nexttile(t); hold(ax3,'on');
plot(weights*100, violationsTrue, 'k-', 'LineWidth',linewidth);
ylabel('Violations (ms)'); xlabel('RT Share (%)');
plot([0 50],[0 0],'r--','LineWidth',linewidth);
axis([0 50 0 10]); box off; ax3.TickDir='out'; ax3.FontSize=fontsize;
yticks(linspace(0, 10, 6))
xticks(linspace(0, 50, 6))
ax3.YColor = "k";
ax3.XColor = "k";
ax3.LineWidth = linewidth;


%% Save figure
outFile = fullfile(pwd, 'Figures', 'figure2.pdf');
exportgraphics(gcf, outFile, 'ContentType', 'vector');


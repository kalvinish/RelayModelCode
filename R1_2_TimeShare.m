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
xx         = linspace(100, 800, 500);  % RT values for CDF evaluation

% Preallocate for speed
modelCDFs = nan(length(xx), numWeights);
gains     = nan(1, numWeights);
violationsTrue = nan(1, numWeights);

%% Load Empirical Quantiles for Plotting
numQuantiles = 10;
empiricalDir = fullfile(pwd, 'EmpiricalData', 'Miller82');
empData = nan(numQuantiles, 3);
for cond = 1:3
    csvFile = fullfile(empiricalDir, sprintf('%d.csv', cond));
    raw = readmatrix(csvFile);
    empData(:, cond) = raw(:,1);
end
probLevels = linspace(0.05, 0.95, numQuantiles);

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

% Unisensory CDFs
uniA = uniCDF(xx, aMU, aLAMBDA);
uniV = uniCDF(xx, vMU, vLAMBDA);
plot(xx, uniA, 'Color','#009FE3', 'LineWidth',1.5);
plot(xx, uniV, 'Color','#3AAA35', 'LineWidth',1.5);

% Relay model curves colored by weight
grey = [228 210 231]/255; black = [137 41 133]/255;
colors = [linspace(grey(1),black(1),numWeights)',...
          linspace(grey(2),black(2),numWeights)',...
          linspace(grey(3),black(3),numWeights)'];
for idx = 1:numWeights
    plot(ax1, xx, modelCDFs(:,idx), 'Color', colors(idx,:), 'LineWidth',1.5);
end

% Reference bounds: Miller and Raab
mxx = linspace(100,800,500);
plot(ax1, mxx, getMillerCDF(mxx, aMU,vMU,aLAMBDA,vLAMBDA), 'r-', 'LineWidth',1.5);
plot(ax1, xx, getRaabCDF(xx, aMU,vMU,aLAMBDA,vLAMBDA), 'r--', 'LineWidth',1.5);

% Optimal-weight relay curve
plot(ax1, xx, getRelayCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA,optW,1-optW,optW,1-optW), 'k-', 'LineWidth',1.5);

% Overlay empirical quantiles
scatter(ax1, empData(:,1), probLevels, 30, 'MarkerEdgeColor','#009FE3', 'MarkerFaceColor','w');
scatter(ax1, empData(:,2), probLevels, 30, 'MarkerEdgeColor','#3AAA35', 'MarkerFaceColor','w');
scatter(ax1, empData(:,3), probLevels, 30, 'MarkerEdgeColor','k',        'MarkerFaceColor','w');

xlabel('Response Time (ms)'); ylabel('Cumulative Probability');
axis([100 700 0 1]); box off;
ax1.TickDir = 'out'; ax1.FontSize = 9;

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
plot(weights*100, gains, 'k-', 'LineWidth',1.5);
ylabel('RSE (ms)');
xlabel('First-stage Weight (%)');
plot([0 weights(end)*100], [getRSE_fromCDF(xx,getRaabCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), getGriceCDF(uniA,uniV))]*[1 1], 'r--','LineWidth',1.5);
axis([0 50 0 80]); box off; ax2.TickDir='out'; ax2.FontSize=9;

% Violation subplot
ax3 = nexttile(t); hold(ax3,'on');
plot(weights*100, violationsTrue, 'k-', 'LineWidth',1.5);
ylabel('Violations (ms)'); xlabel('First-stage Weight (%)');
plot([0 50],[0 0],'r--','LineWidth',1.5);
axis([0 50 0 max(violationsTrue)*1.1]); box off; ax3.TickDir='out'; ax3.FontSize=9;

% Save figure
outFile = fullfile(pwd, 'Figures', 'figure2.pdf');
exportgraphics(gcf, outFile, 'ContentType', 'vector');


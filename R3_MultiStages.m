%% Figure 3: Impact of Number of Relay Stages on CDF, RSE, and Violations
% This script examines how increasing the number of relay (race) stages
% affects the predicted CDFs, redundancy gain (RSE), and Miller-bound violations.

clear; clc; close all;

% Add custom functions to path
addpath(genpath(fullfile(pwd, 'Functions')));

writedata = true;

%% Load Empirical Data and Fitted Parameters from Previous Analysis

paramsFile = fullfile(pwd, 'FittedParameters', 'params_uni_miller_82.mat');
load(paramsFile);

% Extract parameters for convenience
aMU      = params_82(1,1);
vMU      = params_82(2,1);

aLAMBDA  = params_82(1,2);
vLAMBDA  = params_82(2,2);

% Load in empirical data from digitised CDF taken from Figure 1 Miller (1982)
empData_82 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller82', 'miller_82.xlsx'));

%% Define Response-Time Grid and Race Range
xMin   = 0;  xMax   = 700;  nX = 1000;
xx      = linspace(xMin, xMax, nX);
raceMax = 10;  % maximum number of stages to test
raceN   = raceMax;

%% Compute Relay CDFs for 1..raceMax Stages (with RSE and violations)

multiCDF  = nan(nX, raceN);

rseVals   = nan(1,raceN);
violVals  = nan(1,raceN);
rmseVals  = nan(1,raceN);

griceCDF  = getGriceCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
raabCDF   = getRaabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
millerCDF = getMillerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);

parfor k = 1:raceN

    multiCDF(:,k) = getMultiRelayCDF(xx, k, aMU, vMU, aLAMBDA, vLAMBDA);

    % Interpolate predicted CDF at empirical RTs
    predAtEmp = interp1(xx, multiCDF(:,k), empData_82(:,3), 'linear', 'extrap');

    % Redundancy gain (RSE) against race (Grice) model
    rseVals(k)   = getRSE_fromCDF(xx, multiCDF(:,k), griceCDF);
    % Violations against Miller bound
    violVals(k)  = getViolation_fromCDF(xx, multiCDF(:,k), millerCDF);
    % RMSE between empirical and predicted at data points
    diff = empData_82(:,4) - predAtEmp;
    rmseVals(k)  = sqrt(mean(diff.^2));

end


%% Plot Figure 3

f1 = figure();
t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

% Figure 3b - CDF -------------------------------
plotOpts = createPlotOpts([100 700], 3, [0 1], 3);
ax1 = nexttile(t, [1 2]); hold(ax1,'on');

% Plot unisensory CDFs
plot(xx, getUniCDF(xx,aMU,aLAMBDA), 'Color',plotOpts.audCol,'LineWidth',plotOpts.linewidth);
plot(xx, getUniCDF(xx,vMU,vLAMBDA), 'Color',plotOpts.visCol,'LineWidth',plotOpts.linewidth);

% Generate colormap for relay stages (from light grey to black)
lightGrey = [0.8 0.8 0.8];
black    = [0    0    0];
stageColors = [linspace(lightGrey(1),black(1),raceN)', ...
               linspace(lightGrey(2),black(2),raceN)', ...
               linspace(lightGrey(3),black(3),raceN)'];

% Plot each multi-stage CDF
for k = 1:raceN
    plot(xx, multiCDF(:,k), 'Color', stageColors(k,:), 'LineWidth',plotOpts.linewidth);
end

% Plot Race and Miller Bound
plot(xx, raabCDF, 'Color',plotOpts.modelCol,'LineWidth',plotOpts.linewidth, 'LineStyle','--');
plot(xx, millerCDF, 'Color',plotOpts.modelCol,'LineWidth',plotOpts.linewidth, 'LineStyle','-');

% Overlay empirical quantiles
markerOpts = {'Marker', 'o', 'MarkerSize', plotOpts.markersize, 'LineStyle', 'none', 'MarkerFaceColor', plotOpts.markerfacecol, 'LineWidth', plotOpts.markerlinewidth};

plot(ax1, empData_82(:,1), empData_82(:,4), markerOpts{:}, 'Color', plotOpts.audCol, 'DisplayName', 'A');
plot(ax1, empData_82(:,2), empData_82(:,4), markerOpts{:}, 'Color', plotOpts.visCol, 'DisplayName', 'V');
plot(ax1, empData_82(:,3), empData_82(:,4), markerOpts{:}, 'Color', plotOpts.audvisCol, 'DisplayName', 'AV');

xlabel(ax1,'Response Time (ms)'); 
ylabel(ax1,'Cumulative Probability');

% Set axis ticks and other properties
ylim(plotOpts.ylim)
xlim(plotOpts.xlim)
yticks(plotOpts.yticks)
xticks(plotOpts.xticks)

% Figure 3c - RSE -------------------------------
markerOpts = {'Marker', 'o', 'MarkerSize', plotOpts.markersize, 'LineStyle', '-', 'MarkerFaceColor', plotOpts.markerfacecol, 'LineWidth', plotOpts.markerlinewidth};

plotOpts = createPlotOpts([0 10], 6, [0 200], 5);
ax2 = nexttile(t); hold(ax2,'on');
yline(rseVals(1),'--', 'Color', plotOpts.modelCol, 'LineWidth',plotOpts.linewidth);
plot(1, rseVals(1), markerOpts{:},'Color',plotOpts.modelCol);
plot(2:raceN, rseVals(2:end), markerOpts{:},'Color',plotOpts.audvisCol);
xlabel(ax2,'Stages'); ylabel(ax2,'RSE (ms)');
ylim(plotOpts.ylim)
xlim(plotOpts.xlim)
yticks(plotOpts.yticks)
xticks(plotOpts.xticks)

% Figure 3d - Violation -------------------------------
plotOpts = createPlotOpts([0 10], 6, [0 100], 6);
ax3 = nexttile(t); hold(ax3,'on');
yline(violVals(1),'--', 'Color', plotOpts.modelCol, 'LineWidth',plotOpts.linewidth);
plot(1, violVals(1), markerOpts{:},'Color',plotOpts.modelCol);
plot(2:raceN, violVals(2:end), markerOpts{:},'Color',plotOpts.audvisCol);
xlabel(ax3,'Stages'); ylabel(ax3,'Violations (ms)');
ylim(plotOpts.ylim)
xlim(plotOpts.xlim)
yticks(plotOpts.yticks)
xticks(plotOpts.xticks)

allAx = [ax1, ax2, ax3];

set(allAx, ...
    'Box',      'off', ...                
    'TickDir',  plotOpts.tickdir, ...
    'FontSize', plotOpts.fontsize, ...
    'XColor',   plotOpts.axisCol, ...
    'YColor',   plotOpts.axisCol, ...
    'LineWidth',plotOpts.linewidth);

% Save figure
if writedata
    outFile = fullfile(pwd, 'Figures', 'Figure3.pdf');
    exportgraphics(gcf, outFile, 'ContentType', 'vector');
end
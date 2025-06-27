%% Fitting Miller (1982) and Miller (1986) Unisensory Data with Inverse Gaussians
% This script loads response time data, fits inverse Gaussian models to 
% unisensory distributions.

clear; clc; close all;

% Add custom functions to MATLAB path
addpath(genpath(fullfile(pwd, 'Functions')));

writedata = true; % set to true to save data/figures

%% Load Empirical Data

% Load in empirical data from digitised CDF taken from Figure 1 Miller (1982)
empData_82 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller82', 'miller_82.xlsx'));

% Load in empirical data from Table 1 Miller (1986): ID: BD (switch BD for KY for participant KY)
empData_uni_86 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller86', 'miller_86_BD_uni.xlsx'));

%% Fit IG Distributions to Empirical Data

params_82 = nan(3, 2);  % [mu, lambda] for each condition (includes IG fit to AV)
params_86 = nan(2, 2);  % [mu, lambda] for each condition

for cond = 1:3
    % Fit inverse Gaussian to unisensory quantile data
    [mu_opt, lambda_opt] = fitIG_fromCDF(empData_82(:,[cond, 4]));
    params_82(cond,:) = [mu_opt, lambda_opt];
end

% Fit IG distributions using method of moments to unisensory summary statistics (Miller, 1986)
% Auditory params
for cond = 1:2
    [mu_opt, lambda_opt] = fitIG_fromDesc(empData_uni_86(1,cond),...
                                          empData_uni_86(3,cond),...
                                          empData_uni_86(4,cond));
    params_86(cond,:) = [mu_opt, lambda_opt];
end

if writedata
    save(fullfile(pwd, 'FittedParameters', 'params_uni_miller_82.mat'), 'params_82');
    save(fullfile(pwd, 'FittedParameters', 'params_uni_miller_86.mat'), 'params_86');
end

%% Plot CDF for Figure 1b

plotOpts = createPlotOpts([100 700], 3, [0 1], 3);

numQuantiles = size(empData_82, 1);

figure('Visible','on');
hold on;
box off;

% Shaded gain and violation regions
fillArea([empData_82(:,3), getGrice(empData_82(:,1:2))], [0.9 0.9 0.9], 1);
fillArea([empData_82(:,3), getMiller(empData_82(:,1:2))], [0.75 0.75 0.75], 1);

cp = getCP(numQuantiles);

% Plot empirical points and model curves
markerOpts = {'o',...
              'MarkerFaceColor',plotOpts.markerfacecol,...
              'MarkerSize',plotOpts.markersize,...
              'LineStyle','-',...
              'LineWidth',plotOpts.markerlinewidth};

plot(empData_82(:,1), cp, markerOpts{:}, 'Color', plotOpts.audCol, 'DisplayName', 'A');
plot(empData_82(:,2), cp, markerOpts{:}, 'Color', plotOpts.visCol, 'DisplayName', 'V');
plot(empData_82(:,3), cp, markerOpts{:}, 'Color', plotOpts.audvisCol, 'DisplayName', 'AV');
plot(getMiller(empData_82(:,1:2)), cp, '-', 'Color',plotOpts.modelCol, 'LineWidth', plotOpts.linewidth, 'DisplayName', 'RMI');
plot(getRaab(empData_82(:,1:2)),  cp, '--','Color',plotOpts.modelCol, 'LineWidth', plotOpts.linewidth, 'DisplayName', 'Race Model');

xlabel('Response Time (ms)'); ylabel('Cumulative Probability');
legend('Location','Southeast','Box','off');

% Set axis ticks and other properties
xlim(plotOpts.xlim)
ylim(plotOpts.ylim)
yticks(plotOpts.yticks)
xticks(plotOpts.xticks)

ax1 = gca;
ax1.TickDir = plotOpts.tickdir;
ax1.FontSize = plotOpts.fontsize;
ax1.YColor = plotOpts.axisCol;
ax1.XColor = plotOpts.axisCol;
ax1.LineWidth = plotOpts.linewidth;

% Save figure
if writedata
    outFile = fullfile(pwd, 'Figures', 'miller82CDF.pdf');
    exportgraphics(gcf, outFile, 'ContentType', 'vector');
end
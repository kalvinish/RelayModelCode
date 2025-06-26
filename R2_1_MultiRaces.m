%% Figure 3: Impact of Number of Relay Stages on CDF, RSE, and Violations
% This script examines how increasing the number of relay (race) stages
% affects the predicted CDFs, redundancy gain (RSE), and Miller-bound violations.

clear; clc; close all;

% Add custom functions to path
addpath(genpath(fullfile(pwd, 'Functions')));

%% Load Previously Fitted Parameters
% Parameters stored as cell array: {A_dist, V_dist, weightStruct}
paramFile = fullfile(pwd, 'FittedParameters', 'Miller82_Parameters.mat');
load(paramFile, 'params');
A = params{1};  V = params{2};  W = params{3};

% Extract Inverse-Gaussian parameters
aMU     = A.mu;    vMU     = V.mu;
aLAMBDA = A.lambda; vLAMBDA = V.lambda;

%% Define Response-Time Grid and Race Range
xMin   = 0;  xMax   = 700;  nX = 1000;
xx      = linspace(xMin, xMax, nX);
raceMax = 10;  % maximum number of stages to test
raceN   = raceMax;

%% Load Empirical Quantile Data
loadedData = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller82', 'miller.xlsx'));
empData = loadedData(:,1:3);
probLevels = loadedData(:,4);

%% Compute Relay CDFs for 1..raceMax Stages
multiCDF = nan(nX, raceN);
parfor k = 1:raceN
    multiCDF(:,k) = getMultiRelayCDF(xx, k, aMU, vMU, aLAMBDA, vLAMBDA);
end

%% Set up Figure with Tiled Layout
f1 = figure('Units','centimeters','OuterPosition',[0 0 16 5.5]);
t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(t, [1 2]); hold(ax1,'on');

linewidth = 2.5;
fontsize = 14;

% Plot unisensory CDFs
plot(xx, getUniCDF(xx,aMU,aLAMBDA), 'Color','#3AAA35','LineWidth',linewidth);
plot(xx, getUniCDF(xx,vMU,vLAMBDA), 'Color','#009FE3','LineWidth',linewidth);

% Generate colormap for relay stages (from light grey to black)
lightGrey = [0.8 0.8 0.8];
black    = [0    0    0];
stageColors = [linspace(lightGrey(1),black(1),raceN)', ...
               linspace(lightGrey(2),black(2),raceN)', ...
               linspace(lightGrey(3),black(3),raceN)'];

% Plot each multi-stage CDF
for k = 1:raceN
    plot(xx, multiCDF(:,k), 'Color', stageColors(k,:), 'LineWidth',linewidth);
end

% Plot Race and Miller Bound
plot(xx, getRaabCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), 'Color','r','LineWidth',linewidth, 'LineStyle','--');
plot(xx, getMillerCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), 'Color','r','LineWidth',linewidth, 'LineStyle','-');

% Overlay empirical quantiles
scatter(empData(:,1), probLevels, 50, 'MarkerEdgeColor','#3AAA35','MarkerFaceColor','w', 'LineWidth', linewidth);
scatter(empData(:,2), probLevels, 50, 'MarkerEdgeColor','#009FE3','MarkerFaceColor','w', 'LineWidth', linewidth);
scatter(empData(:,3), probLevels, 50, 'MarkerEdgeColor','k','MarkerFaceColor','w', 'LineWidth', linewidth);

xlabel(ax1,'Response Time (ms)'); ylabel(ax1,'Cumulative Probability');
axis(ax1,[100 xMax 0 1]); set(ax1,'Box','off','TickDir','out','FontSize',9,'LineWidth',linewidth);

% Set axis ticks and other properties
yticks(linspace(0, 1, 3))
xticks(linspace(100, 700, 3))

ax1.YColor = "k";
ax1.XColor = "k";
ax1.LineWidth = linewidth;
ax1.FontSize = fontsize;

%% Compute RSE, Violations, and RMSE Across Stages
% Prepare empirical CDF at data points
dataRT      = empData(:,3);            % AV quantile RTs
empiricalF  = probLevels(:);           % CDF levels at those RTs

rseVals     = nan(1,raceN);
violVals    = nan(1,raceN);
rmseVals    = nan(1,raceN);

griceCDF    = getGriceCDF(getUniCDF(xx,aMU,aLAMBDA), getUniCDF(xx,vMU,vLAMBDA));
millerCDF   = getMillerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);

for k = 1:raceN
    % Predicted full-grid CDF
    predCDF = multiCDF(:,k);
    % Interpolate predicted CDF at empirical RTs
    predAtEmp = interp1(xx, predCDF, dataRT, 'linear', 'extrap');

    % Redundancy gain (RSE) against race (Grice) model
    rseVals(k)   = getRSE_fromCDF(xx, predCDF, griceCDF);
    % Violations against Miller bound
    violVals(k)  = getViolation_fromCDF(xx, predCDF, millerCDF);
    % RMSE between empirical and predicted at data points
    diff = empiricalF - predAtEmp;
    rmseVals(k)  = sqrt(mean(diff.^2));
end

%% Plot RSE vs. Number of Stages
ax2 = nexttile(t); hold(ax2,'on');
plot(1:raceN, rseVals, 'o-','Color','k','LineWidth',linewidth,'MarkerFaceColor','w');
yline(rseVals(1),'r--','LineWidth',linewidth);
xlabel(ax2,'Stages'); ylabel(ax2,'RSE (ms)');
set(ax2,'Box','off','TickDir','out','FontSize',fontsize,'LineWidth',linewidth);
xlim([0 raceMax]); ylim([0 200]);

% Set axis ticks and other properties
yticks(linspace(0, 200, 5))
xticks(linspace(0, 10, 6))

ax2.YColor = "k";
ax2.XColor = "k";
ax2.LineWidth = linewidth;
ax2.FontSize = fontsize;

%% Plot Violations vs. Number of Stages
ax3 = nexttile(t); hold(ax3,'on');
plot(1:raceN, violVals, 'o-','Color','k','LineWidth',linewidth,'MarkerFaceColor','w');
yline(violVals(1),'r--','LineWidth',linewidth);
xlabel(ax3,'Stages'); ylabel(ax3,'Violations (ms)');
set(ax3,'Box','off','TickDir','out','FontSize',fontsize,'LineWidth',linewidth);
xlim([0 raceMax]); ylim([0 100]);

% Set axis ticks and other properties
yticks(linspace(0, 100, 6))
xticks(linspace(0, 10, 6))

ax3.YColor = "k";
ax3.XColor = "k";
ax3.LineWidth = linewidth;
ax3.FontSize = fontsize;

%% Save Figure as Vector PDF
outDir = fullfile(pwd,'Figures'); if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(f1, fullfile(outDir,'figure3.pdf'),'ContentType','vector');

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
empDir = fullfile(pwd, 'EmpiricalData', 'Miller82');
numQuantiles = 10;
empData      = nan(numQuantiles, 3);
for cond = 1:3
    file = fullfile(empDir, sprintf('%d.csv', cond));
    tmp  = readmatrix(file);
    empData(:,cond) = tmp(:,1);
end
probLevels = linspace(0.05, 0.95, numQuantiles);

%% Compute Relay CDFs for 1..raceMax Stages
multiCDF = nan(nX, raceN);
parfor k = 1:raceN
    multiCDF(:,k) = getMultiRelayCDF(xx, k, aMU, vMU, aLAMBDA, vLAMBDA);
end

%% Set up Figure with Tiled Layout
f1 = figure('Units','centimeters','OuterPosition',[0 0 16 5.5]);
t = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(t, [1 2]); hold(ax1,'on');

% Plot unisensory CDFs
plot(xx, uniCDF(xx,aMU,aLAMBDA), 'Color','#009FE3','LineWidth',1.5);
plot(xx, uniCDF(xx,vMU,vLAMBDA), 'Color','#3AAA35','LineWidth',1.5);

% Define a colormap using cbrewer2
races = 1:raceN;
colors = cbrewer2('Reds', length(races));

% Generate colormap for relay stages (from light grey to black)
lightGrey = [0.8 0.8 0.8];
black    = [0    0    0];
stageColors = [linspace(lightGrey(1),black(1),raceN)', ...
               linspace(lightGrey(2),black(2),raceN)', ...
               linspace(lightGrey(3),black(3),raceN)'];

% Plot each multi-stage CDF
for k = 1:raceN
    plot(xx, multiCDF(:,k), 'Color', stageColors(k,:), 'LineWidth',1.5);
end

plot(xx, getRaabCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), 'Color','r','LineWidth',1.5, 'LineStyle','--');
plot(xx, getMillerCDF(xx,aMU,vMU,aLAMBDA,vLAMBDA), 'Color','r','LineWidth',1.5, 'LineStyle','-');

% Overlay empirical quantiles
scatter(empData(:,1), probLevels, 30, 'MarkerEdgeColor','#009FE3','MarkerFaceColor','w');
scatter(empData(:,2), probLevels, 30, 'MarkerEdgeColor','#3AAA35','MarkerFaceColor','w');
scatter(empData(:,3), probLevels, 30, 'MarkerEdgeColor','k','MarkerFaceColor','w');

xlabel(ax1,'Response Time (ms)'); ylabel(ax1,'Cumulative Probability');
axis(ax1,[100 xMax 0 1]); set(ax1,'Box','off','TickDir','out','FontSize',9,'LineWidth',1.5);

%% Compute RSE, Violations, and RMSE Across Stages
% Prepare empirical CDF at data points
dataRT      = empData(:,3);            % AV quantile RTs
empiricalF  = probLevels(:);           % CDF levels at those RTs

rseVals     = nan(1,raceN);
violVals    = nan(1,raceN);
rmseVals    = nan(1,raceN);

griceCDF    = getGriceCDF(uniCDF(xx,aMU,aLAMBDA), uniCDF(xx,vMU,vLAMBDA));
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
plot(1:raceN, rseVals, 'o-','Color','k','LineWidth',1.5,'MarkerFaceColor','w');
yline(rseVals(1),'r--','LineWidth',1.5);
xlabel(ax2,'Stages'); ylabel(ax2,'RSE (ms)');
set(ax2,'Box','off','TickDir','out','FontSize',9,'LineWidth',1.5);
xlim([1 raceMax]); ylim([0 max(rseVals)*1.1]);

%% Plot Violations vs. Number of Stages
ax3 = nexttile(t); hold(ax3,'on');
plot(1:raceN, violVals, 'o-','Color','k','LineWidth',1.5,'MarkerFaceColor','w');
yline(violVals(1),'r--','LineWidth',1.5);
xlabel(ax3,'Stages'); ylabel(ax3,'Violations (ms)');
set(ax3,'Box','off','TickDir','out','FontSize',9,'LineWidth',1.5);
xlim([1 raceMax]); ylim([0 max(violVals)*1.1]);

%% Save Figure as Vector PDF
outDir = fullfile(pwd,'Figures'); if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(f1, fullfile(outDir,'figure3.pdf'),'ContentType','vector');

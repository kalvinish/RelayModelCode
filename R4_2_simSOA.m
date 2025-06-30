%% Figure 4: RT-share and Signal Onset Asynchrony
% This script simulates the RSE across a fine grid of SOA values to be plotted
% with the empirical mean RSEs from Miller (1986).

clear; clc; close all;

% Add custom function directory to path
addpath(fullfile(pwd, 'Functions'));

writedata = true; % do you want to write data and figures generated?

%% Load Empirical Data and Fitted Parameters from Previous Analysis

paramsFile = fullfile(pwd, 'FittedParameters', 'params_uni_miller_86.mat');
load(paramsFile); % loads unisensory parameters

bestWfile = fullfile(pwd, 'FittedParameters', 'miller_86_weight.mat');
load(bestWfile); % loads best fitting RT share

% Extract parameters for convenience
params.aMU      = params_86(1,1);
params.vMU      = params_86(2,1);

params.aLAMBDA  = params_86(1,2);
params.vLAMBDA  = params_86(2,2);

params.w = w.opt;

% Load in empirical data from digitised CDF taken from Figure 1 Miller (1982)
empData_86 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller86', 'miller_86_BD_full.xlsx'));

% Row Names
row.soa = 1; row.mean = 2; row.median = 3; row.se = 4; row.rse = 5; row.sd = 6;

% Column Names
col.uniA = 1; col.soa = 2:12; col.uniV = 13;

% Get empirical data
empData.uniA_mean = empData_86(row.mean,col.uniA);
empData.uniV_mean = empData_86(row.mean,col.uniV);
empData.av_mean = empData_86(row.mean, col.soa);
empData.rse = empData_86(row.rse, col.soa);
empData.lags = empData_86(row.soa, col.soa);
empData.nTrials = 400;

%% Get time for component stages

realA.firstStageMean = params.aMU * params.w;
realA.firstStageSD = sqrt((realA.firstStageMean^3)/(params.aLAMBDA * params.w));
realA.secondStageMean = params.aMU * (1-params.w);
realA.secondStageSD = sqrt((realA.secondStageMean^3)/(params.aLAMBDA * (1-params.w)));

realV.firstStageMean = params.vMU * params.w;
realV.firstStageSD = sqrt((realV.firstStageMean^3)/(params.vLAMBDA * params.w));
realV.secondStageMean = params.vMU * (1-params.w);
realV.secondStageSD = sqrt((realV.secondStageMean^3)/(params.vLAMBDA * (1-params.w)));

%% 4. SIMULATE SOA EXPERIMENT
% Monte Carlo simulation over fine SOA grid
repN   = 100000; % Number of repetitions

% Finer SOA grid for simulation (ms)
simLags   = linspace(-167, 167, 100);

% Run SOA simulation
[meanA, meanV, meanRelay, meanRace, rseRelay, rseRace] = runSOAsim(empData, params, simLags, repN);

% Get 95% CI
[relayCI_low, relayRSE_mean, relayCI_high] = computeCI(rseRelay, 2);
[raceCI_low, raceRSE_mean,  raceCI_high]  = computeCI(rseRace,  2);
[relayMeanCI_low, relayMean_mean, relayMeanCI_high]  = computeCI(meanRelay,  2);
[raceMeanCI_low, raceMean_mean, raceMeanCI_high]  = computeCI(meanRace,  2);

%% Plot Results

plotOpts = createPlotOpts([-200 200], 5, [200 400], 5);

lagsA = zeros(length(simLags), 1);
lagsA(simLags > 0) = simLags(simLags > 0);
uniA_plusLag = (empData.uniA_mean + lagsA);

lagsV = zeros(length(simLags), 1);
lagsV(simLags < 0) = abs(simLags(simLags < 0));
uniV_plusLag = (empData.uniV_mean + lagsV);

figure();
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel 1: Mean RT vs. SOA
ax1 = nexttile(t);
hold on
% CI shading for relay model medians
xShade = [simLags, fliplr(simLags)];
yShade = [relayMeanCI_low', fliplr(relayMeanCI_high')];
fill(xShade, yShade, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;
% Plot unisensory and multisensory curves
plot(simLags, uniA_plusLag, '--', 'LineWidth', plotOpts.linewidth, 'DisplayName', 'A', 'Color',plotOpts.audCol);
plot(simLags, uniV_plusLag, '--', 'LineWidth', plotOpts.linewidth, 'DisplayName', 'V', 'Color',plotOpts.visCol);
plot(simLags, relayMean_mean, '-', 'LineWidth', plotOpts.linewidth, 'DisplayName', 'Relay', 'Color',plotOpts.audvisCol);
plot(simLags, raceMean_mean, '--','LineWidth', plotOpts.linewidth, 'DisplayName', 'Race', 'Color',plotOpts.modelCol);
% Overlay empirical mean RTs
markerOpts = {'Marker', 'o', 'MarkerSize', plotOpts.markersize, 'LineStyle', 'none', 'MarkerFaceColor', plotOpts.markerfacecol, 'LineWidth', plotOpts.markerlinewidth};
plot(empData.lags, empData.av_mean, markerOpts{:}, 'Color', plotOpts.audvisCol, 'DisplayName', 'Miller (1986)');
xlabel('SOA (ms)'); ylabel('Response Time (ms)');
ylim(plotOpts.ylim)
xlim(plotOpts.xlim)
yticks(plotOpts.yticks);
xticks(plotOpts.xticks);

% Panel 2: RSE vs. SOA
ax2 = nexttile(t);
hold on;
% CI shading for relay and race models
fill(xShade, [relayCI_low', fliplr(relayCI_high')], 'k', 'FaceAlpha',0.2, 'EdgeColor','none');
fill(xShade, [raceCI_low', fliplr(raceCI_high')], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
% Plot mean RSE curves
plot(simLags, relayRSE_mean, 'LineWidth',plotOpts.linewidth,'DisplayName','Relay', 'Color',plotOpts.audvisCol);
plot(simLags, raceRSE_mean, '--', 'LineWidth',plotOpts.linewidth,'DisplayName','Race', 'Color', plotOpts.modelCol);
% Overlay empirical RSEs
plot(empData.lags, empData.rse, markerOpts{:}, 'Color', plotOpts.audvisCol, 'DisplayName', 'Miller (1986)');
yline(0,'--k'); 
xlabel('SOA (ms)'); ylabel('RSE (ms)');
plotOpts = createPlotOpts([-200 200], 5, [-20 80], 6);
ylim(plotOpts.ylim)
xlim(plotOpts.xlim)
yticks(plotOpts.yticks);
xticks(plotOpts.xticks);

allAx = [ax1, ax2];
set(allAx, ...
    'Box',      'off', ...                
    'TickDir',  plotOpts.tickdir, ...
    'FontSize', plotOpts.fontsize, ...
    'XColor',   plotOpts.axisCol, ...
    'YColor',   plotOpts.axisCol, ...
    'LineWidth',plotOpts.linewidth);

%%
% Save figure
if writedata
    outFile = fullfile(pwd, 'Figures', 'Figure4.pdf');
    exportgraphics(gcf, outFile, 'ContentType', 'vector');
end

%% Local Functions

function [meanA, meanV, meanRelay, meanRace, rseRelay, rseRace] = runSOAsim(empData, params, simLags, repN)
    % Ensure simLags is a row vector for parfor slicing
    simLags = simLags(:)';
    nSimLags = numel(simLags);

    % Preallocate output arrays
    meanA      = zeros(nSimLags, repN);
    meanV      = zeros(nSimLags, repN);
    meanRelay  = zeros(nSimLags, repN);
    meanRace   = zeros(nSimLags, repN);
    rseRelay   = zeros(nSimLags, repN);
    rseRace    = zeros(nSimLags, repN);

    % Create and display the waitbar
    h = waitbar(0, 'Simulation Progress');

    % DataQueue for progress updates
    dq = parallel.pool.DataQueue;
    numCompleted = 0;
    afterEach(dq, @progressUpdate);

    parfor rep = 1:repN
        % Local temp storage to avoid slicing issues in nested loops
        localMeanA     = zeros(nSimLags,1);
        localMeanV     = zeros(nSimLags,1);
        localMeanRelay = zeros(nSimLags,1);
        localMeanRace  = zeros(nSimLags,1);
        localRseRelay  = zeros(nSimLags,1);
        localRseRace   = zeros(nSimLags,1);

        for iLag = 1:nSimLags % loop over lags
            lag = simLags(iLag);
            lagA = max(0, lag); % set the auditory lag
            lagV = max(0, -lag); % set the visual lag

            rtA = empData.uniA_mean + lagA; % set the unisensory A RT
            rtV = empData.uniV_mean + lagV; % set the unisensory V RT

            % Simulate Relay model with best fitting RT share
            rtRelay = getRelayLagRND(empData.nTrials, params.aMU, params.vMU, ...
                                     params.aLAMBDA, params.vLAMBDA, ...
                                     params.w, 1-params.w, params.w, 1-params.w, ...
                                     lag, 1);

            % Simulate Race model where RT share = 0
            rtRaceModel = getRelayLagRND(empData.nTrials, params.aMU, params.vMU, ...
                                         params.aLAMBDA, params.vLAMBDA, ...
                                         0, 1, 0, 1, ...
                                         lag, 1);

            uniMin = min(rtA, rtV); % calculate the 'winning' unisensory

            localMeanA(iLag)     = mean(rtA);
            localMeanV(iLag)     = mean(rtV);
            localMeanRelay(iLag) = mean(rtRelay);
            localMeanRace(iLag)  = mean(rtRaceModel);
            localRseRelay(iLag)  = uniMin - mean(rtRelay); % RSE = fastest uni - multi
            localRseRace(iLag)   = uniMin - mean(rtRaceModel);
        end
        % Assign to sliced outputs
        meanA(:,rep)     = localMeanA;
        meanV(:,rep)     = localMeanV;
        meanRelay(:,rep) = localMeanRelay;
        meanRace(:,rep)  = localMeanRace;
        rseRelay(:,rep)  = localRseRelay;
        rseRace(:,rep)   = localRseRace;

        send(dq, rep);
    end

    close(h);

    function progressUpdate(~)
        numCompleted = numCompleted + 1;
        waitbar(numCompleted/repN, h, sprintf('Completed %d of %d reps', numCompleted, repN));
    end
end

function [qLow, meanVal, qHigh] = computeCI(data, dim)
    % Get quantile confidence intervals (95%)
    qLow  = quantile(data,0.025,dim);
    qHigh = quantile(data,0.975,dim);
    meanVal  = mean(data,dim);
end
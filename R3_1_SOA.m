%% SIMULATE SOA DATA USING EMPIRICAL DATA FROM MILLER (1986)
% This script simulates stimulus onset asynchrony (SOA) reaction time (RT)
% data using empirical measurements from Miller (1986), fits optimal
% weight parameters by minimising RMSE, and visualises the results over both
% empirical and finer simulation grids.

%% 0. INITIAL SETUP
clearvars;    % Remove all variables from workspace
clc;          % Clear command window
close all;    % Close all figure windows
addpath(fullfile(pwd, 'Functions'));  % Add custom functions to path

%% 1. FIT INVERSE GAUSSIAN PARAMETERS
% Fit auditory (A) and visual (V) RT distributions using descriptive stats
% Inputs: mean (ms), median, SE, number of samples
[muA, lambdaA] = fitIG_fromDesc(231, 219, 2.8, 400);  % Auditory params
[muV, lambdaV] = fitIG_fromDesc(348, 326, 4.6, 400);  % Visual params

mean_rtA = 231;
mean_rtV = 348;

% Uncomment below for participant KY data
% [muA, lambdaA, res1] = fitIGfromDesc(0.211, 0.003, 0.193, 400);
% [muV, lambdaV, res2] = fitIGfromDesc(0.282, 0.0031, 0.266, 400);

%% 2. LOAD EMPIRICAL SOA DATA
% Empirical SOA lags (ms) from Miller (1986)
empiricalLags   = [-167, -133, -100, -67, -33, 0, 33, 67, 100, 133, 167];
empiricalMeans  = [234, 230, 227, 228, 221, 217, 238, 263, 277, 298, 316];
empiricalMedians= [223, 222, 221, 223, 219, 215, 233, 259, 274, 296, 316];
empiricalSE     = [2.9, 2.0, 2.0, 1.6, 1.4, 1.4, 1.4, 1.3, 1.5, 1.6, 1.7];
empiricalRSE    = [-3, 1, 4, 3, 10, 14, 26, 35, 54, 50, 32];

% Uncomment below for participant KY data
% millerMeans = [216 217 214 218 215 208 237 249 256 273 278];
% millerMedians = [196 194 192 192 198 190 216 235 244 264 268];
% millerMeanSE = [3.7 3.8 3.9 3.8 3.3 3.2 3.1 2.9 2.3 2.7 3.1];
% millerRSE = [-5 -6 -3 -7 -4 3 7 29 26 9 4];

% Finer SOA grid for simulation (ms)
simLags   = linspace(-167, 167, 200);
nEmpLags  = numel(empiricalLags);
nSimLags  = numel(simLags);

% Number of trials per simulation
nTrials   = 400;

%% 3. OPTIMISE WEIGHT USING RMSE
% Optimise weight w in [0,0.5] to minimise RMSE between simulated and
% empirical RSE values over empirical lags.
nOptimSamples = 100000;  % Samples per RMSE evaluation
nRepeats      = 10;      % Number of optimisation repeats
wEstimates    = zeros(nRepeats,1);
fvals         = zeros(nRepeats,1);

parfor rep = 1:nRepeats
    objFun = @(w) computeRMSE(w, nOptimSamples, empiricalLags, muA, muV, lambdaA, lambdaV, empiricalRSE, mean_rtA, mean_rtV);
    [wEstimates(rep), fvals(rep)] = fminbnd(objFun, 0, 0.5);
    disp(rep)
end

% Aggregate results
wOpt = mean(wEstimates);
wStd = std(wEstimates);
fprintf('Optimal weight: %.4f ± %.4f (mean ± std)\n', wOpt, wStd);

% Define weights for relay model stages based on wOpt
w1 = [wOpt, 1 - wOpt];
w2 = w1;

%% 4. SIMULATE SOA EXPERIMENT
% Monte Carlo simulation over fine SOA grid
repN   = 100;             % Number of repetitions
meanA      = nan(nSimLags, repN);
meanV      = nan(nSimLags, repN);
meanRelay  = nan(nSimLags, repN);
meanRace   = nan(nSimLags, repN);
rseRelay   = nan(nSimLags, repN);
rseRace    = nan(nSimLags, repN);

% Progress indicator
hWait = waitbar(0, 'Simulating SOA trials...');
for rep = 1:repN
    for iLag = 1:nSimLags
        lag = simLags(iLag);
        % Determine modality-specific delays (ms)
        lagA = max(0,  lag);
        lagV = max(0, -lag);
        
        % Generate unisensory RTs with delay
        % rtA = getUniRND(nTrials, muA, lambdaA, 1) + lagA;
        % rtV = getUniRND(nTrials, muV, lambdaV, 1) + lagV;

        rtA = mean_rtA + lagA;
        rtV = mean_rtV + lagV;
        
        % Relay model RTs (weighted integration)
        rtRelay = getRelayLagRND(nTrials, muA, muV, lambdaA, lambdaV, w1(1), w1(2), w2(1), w2(2), lag, 1);
        % Race model RTs (independent racing)
        rtRaceModel = getRelayLagRND(nTrials, muA, muV, lambdaA, lambdaV, 0, 1, 0, 1, lag, 1);
        
        % Compute minimal unisensory mean RT
        uniMin = min(mean(rtA), mean(rtV));
        
        % Store means and redundancy gains
        meanA(iLag,rep)      = mean(rtA);
        meanV(iLag,rep)      = mean(rtV);
        meanRelay(iLag,rep)  = mean(rtRelay);
        meanRace(iLag,rep)   = mean(rtRaceModel);
        rseRelay(iLag,rep)   = uniMin - mean(rtRelay);
        rseRace(iLag,rep)    = uniMin - mean(rtRaceModel);
    end
    waitbar(rep/repN, hWait);
end
close(hWait);

%% 5. COMPUTE SUMMARY STATISTICS
% Compute 95% CI for RSEs and Means
[relayCI_low, relayCI_high] = computeQuantileCI(rseRelay, 2);
[raceCI_low,  raceCI_high]  = computeQuantileCI(rseRace,  2);
[relayMeanCI_low, relayMeanCI_high]  = computeQuantileCI(meanRelay,  2);

% Compute means across repetitions
meanA_all     = mean(meanA, 2);
meanV_all     = mean(meanV, 2);
meanRelay_all = mean(meanRelay, 2);
meanRace_all  = mean(meanRace, 2);
rseRelay_all  = mean(rseRelay, 2);
rseRace_all   = mean(rseRace, 2);

%%
test_rse = min(meanA_all, meanV_all) - meanRelay_all

%% 6. PLOT RESULTS
figure();
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

aud_col = '#3AAA35';
vis_col = '#009FE3';
relay_col = '#706F6F';

lagsA = zeros(length(simLags), 1);
lagsA(simLags > 0) = simLags(simLags > 0);
aData = (mean_rtA + lagsA);

lagsV = zeros(length(simLags), 1);
lagsV(simLags < 0) = abs(simLags(simLags < 0));
vData = (mean_rtV + lagsV);

% Panel 1: Mean RT vs. SOA
ax1 = nexttile(t);
% CI shading for relay model medians
xShade = [simLags, fliplr(simLags)];
yShade = [relayMeanCI_low', fliplr(relayMeanCI_high')];
fill(xShade, yShade, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;
% Plot unisensory and multisensory curves
plot(simLags, aData, '--', 'LineWidth', 1.5, 'DisplayName', 'A', 'Color',aud_col);
plot(simLags, vData, '--b', 'LineWidth', 1.5, 'DisplayName', 'V', 'Color',vis_col);
plot(simLags, meanRelay_all, '-', 'LineWidth', 1.5, 'DisplayName', 'Relay', 'Color',relay_col);
plot(simLags, meanRace_all,   '--','LineWidth', 1.5, 'DisplayName', 'Race', 'Color', 'r');
% Overlay empirical mean RTs
scatter(empiricalLags, empiricalMeans, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'DisplayName', 'Miller (1986)')
xlabel('SOA (ms)'); ylabel('RT (ms)');
ylim([200 400])
box off;

% Panel 2: RSE vs. SOA
ax2 = nexttile(t);
% CI shading for relay and race models
fill(xShade, [relayCI_low', fliplr(relayCI_high')], [0.4 0.4 0.4], 'FaceAlpha',0.2, 'EdgeColor','none');
hold on;
fill(xShade, [raceCI_low', fliplr(raceCI_high')], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
% Plot mean RSE curves
plot(simLags, rseRelay_all, 'LineWidth',1.5,'DisplayName','Relay', 'Color',relay_col);
plot(simLags, rseRace_all, '--', 'LineWidth',1.5,'DisplayName','Race', 'Color', 'r');
% Overlay empirical RSEs
scatter(empiricalLags, empiricalRSE,50,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'DisplayName', 'Miller (1986)')
yline(0,'--k'); 
xlabel('SOA (ms)'); ylabel('RSE (ms)');
box off;

%% 7. SAVE FIGURE
exportgraphics(gcf, fullfile(pwd, 'Figures', 'figure4.pdf'), 'ContentType','vector');

%% LOCAL FUNCTIONS

function rmse = computeRMSE(w, nSamples, lags, muA, muV, lamA, lamV, empiricalRSE, mean_rtA, mean_rtV)
    if numel(w)>1, w1=w; w2=w; else, w1=[w,1-w]; w2=w1; end
    RSEsim = arrayfun(@(lag) simulateRSE(lag, nSamples, muA, muV, lamA, lamV, w1, w2, mean_rtA, mean_rtV), lags);
    rmse = sqrt(mean((empiricalRSE - RSEsim).^2));
end

function rse = simulateRSE(lag, n, muA, muV, lamA, lamV, w1, w2, mean_rtA, mean_rtV)
    lagA = max(0, lag); lagV = max(0, -lag);
    % rtA = getUniRND(n, muA, lamA, 1) + lagA;
    % rtV = getUniRND(n, muV, lamV, 1) + lagV;
    rtA = mean_rtA + lagA;
    rtV = mean_rtV + lagV;

    rtRelay = getRelayLagRND(n, muA, muV, lamA, lamV, w1(1),w1(2),w2(1),w2(2),lag,1);
    rse = min(mean(rtA),mean(rtV)) - mean(rtRelay);
end

function [qLow, qHigh] = computeQuantileCI(data, dim)
    qLow  = quantile(data,0.025,dim);
    qHigh = quantile(data,0.975,dim);
end

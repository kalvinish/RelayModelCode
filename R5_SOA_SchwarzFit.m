%% Figure 2: RT-share Manipulation and Violation Effects
% This script visualises how varying the first-stage RT-share in the relay model
% affects predicted CDFs, RSE, and Miller bound violations.

clear; clc; close all;

% Add custom function directory to path
addpath(fullfile(pwd, 'Functions'));

writedata = true;

%% Load Empirical Data and Fitted Parameters from Previous Analysis

% Load in empirical data from digitised CDF taken from Figure 1 Miller (1982)
empData_86 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller86', 'miller_86_BD_full.xlsx'));

soaRow = 1; % multi - row holding soa values
meanRow = 2; % multi - row holding means
medianRow = 3; % multi - row holding medians
seRow = 4; % multi - row holding SE
sdRow = 6;

% Change SOA for unisensory conditions to be -Inf and Inf
empData_86(empData_86(:,:) == -10000) = -Inf;
empData_86(empData_86(:,:) == 10000) = Inf;

%% Fitting

lags = empData_86(soaRow,:);
empN = 400;
empMean = empData_86(meanRow,:);
empMedian = empData_86(medianRow,:);
empSE = empData_86(seRow,:);
empSD = empSE .* sqrt(empN);

nBoot = 10;
simN = 400;

% Parameter Bounds
uBound = [Inf Inf Inf Inf 0.5 0.5];
lBound = [eps eps eps eps 0 0];

% Start Values
[start_muA, start_lambdaA] = fitIG_fromDesc(empMean(1), empMedian(1), empSE(1), empN);
[start_muV, start_lambdaV] = fitIG_fromDesc(empMean(end), empMedian(end), empSE(end), empN);

startVals = [start_muA, start_muV, start_lambdaA, start_lambdaV, 0.25, 0.2];

% Objective function handle
objFun = @(p) computeSSD(p, simN, nBoot, lags, empMean, empSD);

opts = optimoptions('fmincon', ...
  'Display',                'iter', ...
  'Algorithm',              'interior-point', ...
  'MaxIterations',          5000, ...
  'MaxFunctionEvaluations', 1e5, ...
  'OptimalityTolerance',    1e-10, ...
  'StepTolerance',          1e-10);

problem = createOptimProblem('fmincon', ...
    'x0',        startVals, ...
    'objective', @(p) objFun(p), ...
    'lb',        lBound, ...
    'ub',        uBound, ...
    'options',   opts);

% Use GlobalSearch
gs = GlobalSearch();
[estP,fval] = run(gs, problem);

%% Simulate Best Fitting Data

% Number of bootstrap replicates
nBoot = 500;  

% Pre‐allocate
predMeanBoot = nan(nBoot, numel(lags));
predSDBoot   = nan(nBoot, numel(lags));

% Bootstrap predicted summaries
parfor b = 1:nBoot
    [m,s] = getSummaryRelay(estP, lags, empN);
    predMeanBoot(b,:) = m;
    predSDBoot(b,:)   = s;
end

% Compute means and 95% quantile‐CIs
mean_predMean = mean(predMeanBoot,1);
ci_predMean   = quantile(predMeanBoot, [0.025 0.975]);

mean_predSD   = mean(predSDBoot,1);
ci_predSD     = quantile(predSDBoot,   [0.025 0.975]);

%% Plot Fits
% Convert lags to percentages if needed
x = 1:length(lags);

% Create figure
figure('Name','Mean & SD with 95% CI','Color','w');

% --- Mean subplot ---
subplot(2,1,1);
hold on; box off;

% empirical
plot(x, empMean, 'ko-', 'LineWidth',1.5, 'DisplayName','Empirical Mean');

% predicted mean + CI (errorbars)
errorbar( x, mean_predMean, ...
          mean_predMean - ci_predMean(1,:), ...
          ci_predMean(2,:) - mean_predMean, ...
          'rs--','LineWidth',1.5,'MarkerSize',6,...
          'DisplayName','Predicted Mean' );

xlabel('SOA (ms)');
ylabel('Mean RT (ms)');
legend('Location','best');
title('Mean RT: Empirical vs. Predicted');

% --- SD subplot ---
subplot(2,1,2);
hold on; box off;

plot(x, empSD, 'ko-', 'LineWidth',1.5, 'DisplayName','Empirical SD');

errorbar( x, mean_predSD, ...
          mean_predSD - ci_predSD(1,:), ...
          ci_predSD(2,:) - mean_predSD, ...
          'rs--','LineWidth',1.5,'MarkerSize',6,...
          'DisplayName','Predicted SD' );

xlabel('SOA (ms)');
ylabel('SD RT (ms)');
legend('Location','best');
title('RT SD: Empirical vs. Predicted');

%% Functions

function SSD = computeSSD(params, simN, nBoot, lags, empMean, empSD)
    
    parfor rep = 1:nBoot
        
        [pMean, pSD] = getSummaryRelay(params, lags, simN);
        predMean(rep,:) = pMean;
        predSD(rep,:) = pSD;
    end

    dMean = mean(predMean) - empMean;
    dSD = mean(predSD) - empSD;

    % SSD = sum([dMean.^2, dSD.^2]);
    SSD = sum(dSD.^2);

end

function [predMean, predSD] = getSummaryRelay(params, lags, simN)

    muA  = params(1);
    muV  = params(2);
    lamA = params(3);
    lamV = params(4);
    w1   = params(5);
    w2   = params(6);

    predMean = nan(1, numel(lags));
    predSD  = nan(1, numel(lags));

    for i = 1:numel(lags)
        rtRelay = getRelayLagRND(simN, muA, muV, lamA, lamV, ...
                                 w1, 1-w1, w2, 1-w2, lags(i), 1);
        predMean(i) = mean(rtRelay);
        predSD(i)   = std(rtRelay);
    end

end
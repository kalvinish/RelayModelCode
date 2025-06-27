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

%% Fitting with Bayesian Optimization

simN = 1000000;

%%— Define your threshold (e.g. we want SSD < 5) —%%
ssdThreshold = 5;

% Wrap computeSSD into a function that takes a table of optimVars
funBO = @(optVars) computeSSD(...
    [optVars.muA,  optVars.muV,  optVars.lamA,  optVars.lamV, ...
     optVars.w1,   optVars.w2], ...
    simN, lags, empMean, empSD);

% Define optimization variables and their bounds
vars = [
    optimizableVariable('muA',   [eps, 10000])
    optimizableVariable('muV',   [eps, 10000])
    optimizableVariable('lamA',  [eps, 100000])
    optimizableVariable('lamV',  [eps, 100000])
    optimizableVariable('w1',    [0, 0.5])
    optimizableVariable('w2',    [0, 0.5])
];

%%— Updated bayesopt call —%%
results = bayesopt(funBO, vars, ...
    'IsObjectiveDeterministic', false, ...
    'UseParallel',               true, ...
    'NumSeedPoints',            30, ...                   % more initial coverage
    'MaxObjectiveEvaluations',  200, ...                  % more total evaluations
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'OutputFcn',               {@stopIfReached}, ...      % our early-stop hook
    'Verbose',                  1, ...
    'PlotFcn',                 {@plotObjectiveModel, @plotAcquisitionFunction} ...
);


% Extract best parameters
best = results.XAtMinObjective;
estP_bo = [best.muA, best.muV, best.lamA, best.lamV, best.w1, best.w2];
fval_bo = results.MinObjective;

fprintf('\nBayesOpt best SSD = %g\n', fval_bo);
disp('Best parameters:');
disp(estP_bo);


%% Simulate Best Fitting Data

% Number of bootstrap replicates
nBoot = 500;  

% Pre‐allocate
predMeanBoot = nan(nBoot, numel(lags));
predSDBoot   = nan(nBoot, numel(lags));

% Bootstrap predicted summaries
parfor b = 1:nBoot
    [m,s] = getSummaryRelay(estP_bo, lags, empN);
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

%%— Make an OutputFcn to halt when MinObjective ≤ threshold —%%
function stop = stopIfReached(results, state)
    stop = false;
    if strcmp(state, 'iteration')
        if results.MinObjective <= 5
            fprintf('>>> Early stopping: reached SSD = %.4f\n', results.MinObjective);
            stop = true;
        end
    end
end

function SSD = computeSSD(params, simN, lags, empMean, empSD)
   
    [pMean, pSD] = getSummaryRelay(params, lags, simN);

    dMean = pMean - empMean;
    dSD = pSD - empSD;

    % SSD = sqrt(mean([dMean.^2, dSD.^2]));
    SSD = sqrt(mean(dMean.^2));

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
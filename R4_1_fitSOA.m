%% Figure 4: RT-share and Signal Onset Asynchrony
% This script fits the RT share parameter to data from Miller (1986) which
% includes an SOA manipulation. Data is fitted using bayesopt in the paper,
% but here we make it an option to fit using fminbnd.

clear; clc; close all;

% Add custom function directory to path
addpath(fullfile(pwd, 'Functions'));

optim_mode = 'bayesopt'; % Options: {'bayesopt', 'fminbnd'}
writedata = true;

%% Load Empirical Data and Fitted Parameters from Previous Analysis

paramsFile = fullfile(pwd, 'FittedParameters', 'params_uni_miller_86.mat');
load(paramsFile);

% Extract parameters for convenience
params.aMU      = params_86(1,1);
params.vMU      = params_86(2,1);

params.aLAMBDA  = params_86(1,2);
params.vLAMBDA  = params_86(2,2);

% Load in empirical data from digitised CDF taken from Figure 1 Miller (1982)
empData_86 = readmatrix(fullfile(pwd, 'EmpiricalData', 'Miller86', 'miller_86_BD_full.xlsx'));

% Row Names
row.soa = 1; row.mean = 2; row.median = 3; row.se = 4; row.rse = 5; row.sd = 6;

% Column Names
col.uniA = 1; col.soa = 2:12; col.uniV = 13;

% Get empirical means from data
empData.uniA_mean = empData_86(row.mean,col.uniA);
empData.uniV_mean = empData_86(row.mean,col.uniV);
empData.av_mean = empData_86(row.mean, col.soa);
empData.rse = empData_86(row.rse, col.soa);
empData.lags = empData_86(row.soa, col.soa);
empData.nTrials = 400;

nSamples = 100000;  % Monte Carlo draws per simulation

%% Optimise RT share using bayesopt or fminbnd

switch lower(optim_mode)
    case 'bayesopt'
        R        = 100;      % Number of replicates per Bayesian step
        baseSeed = 42;      % Base seed for CRN across replicates
        vars = optimizableVariable('w',[0,0.5],'Type','real');
        objFcn = @(x) replicatedRMSE(x, params, empData, nSamples, R, baseSeed);
        
        results = bayesopt(objFcn, vars, ...
            'MaxObjectiveEvaluations', 100,...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'IsObjectiveDeterministic', false, ...
            'UseParallel', true, ...
            'Verbose', 1,...
            'PlotFcn','all');
        
        % Best result
        w.opt = results.XAtMinObjective.w;
        w.rmse = results.MinObjective;
        fprintf('Best w = %.4f  (RMSE ≈ %.4f)\n', w.opt, w.rmse);
    
        if writedata
            outFile = fullfile(pwd, 'Figures', 'weightFit_miller86.pdf');
            fig1 = figure(1);
            exportgraphics(fig1, outFile, 'ContentType', 'vector');
        end

    case 'fminbnd'
        % Optimise weight RT share in [0,0.5] to minimise RMSE between simulated and
        % empirical RSE values over empirical lags.
        nRepeats      = 10;      % Number of optimisation repeats
        
        [wEstimates, fvals] = optimiseRTshare(params, empData, nSamples, nRepeats);
        
        % Aggregate results
        w.opt = mean(wEstimates);
        w.opt_sd = std(wEstimates);
        w.rmse = mean(fvals);
        w.rmse_sd = std(fvals);
        fprintf('Optimal weight: %.4f ± %.4f (mean ± std)\nRMSE: %.4f ± %.4f ms\n', w.opt, w.opt_sd, w.rmse, w.rmse_sd);

end

%% Save best fitting RT share

if writedata
    outFile = fullfile(pwd, 'FittedParameters', 'miller_86_weight');
    save(outFile, 'w');
end

%% Functions

function obj = replicatedRMSE(x, params, empData, nSamples, R, baseSeed)
    sims = zeros(R,1);
    for k = 1:R
        % Set seed for common random numbers: same stream index across x
        rng(baseSeed + k, 'twister');
        sims(k) = computeRMSE(x.w, params, empData, nSamples);
    end
    obj = mean(sims);
end

% Add a second function for fminbnd estimation with progress
function [wEstimates, fvals] = optimiseRTshare(params, empData, nOptimSamples, nRepeats)
    % Preallocate outputs
    wEstimates = zeros(nRepeats,1);
    fvals      = zeros(nRepeats,1);

    % Create and display the waitbar
    h = waitbar(0, 'Estimation Progress');

    % DataQueue for progress updates
    dq = parallel.pool.DataQueue;
    numCompleted = 0;
    afterEach(dq, @progressUpdate);

    parfor rep = 1:nRepeats
        % Define objective
        objFun = @(w) computeRMSE(w, params, empData, nOptimSamples);
        % Run bounded minimization
        [wEstimates(rep), fvals(rep)] = fminbnd(objFun, 0, 0.5);
        % Signal completion
        send(dq, rep);
    end

    close(h);

    function progressUpdate(~)
        numCompleted = numCompleted + 1;
        waitbar(numCompleted/nRepeats, h, sprintf('Completed %d of %d repeats', numCompleted, nRepeats));
    end
end

function rmse = computeRMSE(w, params, empData, nSamples)   
    RSEsim = arrayfun(@(lag) simulateRSE(w, params, empData, lag, nSamples), empData.lags);
    rmse = sqrt(mean((empData.rse - RSEsim).^2));
end

function rse = simulateRSE(w, params, empData, lag, nSamples)

    lagA = max(0, lag); lagV = max(0, -lag);

    rtA = empData.uniA_mean + lagA;
    rtV = empData.uniV_mean + lagV;

    rtRelay = getRelayLagRND(nSamples,...
                             params.aMU, params.vMU,...
                             params.aLAMBDA, params.vLAMBDA,...
                             w, 1-w, w, 1-w,...
                             lag,1);
    rse = min(rtA,rtV) - mean(rtRelay);
end
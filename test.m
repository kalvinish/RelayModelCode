%% parameterRecovery_v3.m
% Fast, reliable parameter‐recovery for the relay model

clear; clc; close all;
addpath(fullfile(pwd,'Functions'));

%% 1) Design & settings
lags   = [-Inf, -50, 0, 50, Inf];
empN   = 400;     % # samples per “empirical” dataset
nRecov = 1;      % # recovery reps
simN   = 1e5;     % draws per lag when simulating summary (large → low noise)

trueP  = [300, 320, 7000, 7500, 0.2];   % ground truth

rt_a = getUniRND(empN, trueP(1), trueP(3));
rt_v = getUniRND(empN, trueP(2), trueP(4));

%% 2) Fit bounds & base start
lb = [eps, eps,   eps,   eps,   0];
ub = [Inf, Inf, Inf, Inf, 0.5];
[muA0, lamA0] = fitIG_fromDesc(mean(rt_a),median(rt_a),std(rt_a), empN);
[muV0, lamV0] = fitIG_fromDesc(mean(rt_v),median(rt_v),std(rt_v), empN);
x0 = [muA0, muV0, lamA0, lamV0, 0.25];

%% 3) fmincon options

% wrap SSD with fixed simN & lags
ssdFun = @(p, m, s) computeSSD(p, lags, m, s, simN);

%% 4) Recovery
recovered = nan(nRecov,5);
fvals      = nan(nRecov,1);

[synMean, synSD] = getSummaryRelay(trueP, lags, simN);

% 1) tighter fmincon options
opts = optimoptions('fmincon', ...
  'Display',                'iter', ...
  'Algorithm',              'interior-point', ...
  'MaxIterations',          5000, ...
  'MaxFunctionEvaluations', 1e5, ...
  'OptimalityTolerance',    1e-10, ...
  'StepTolerance',          1e-10);

problem = createOptimProblem('fmincon', ...
    'x0',        x0, ...
    'objective', @(p) ssdFun(p,synMean,synSD), ...
    'lb',        lb, ...
    'ub',        ub, ...
    'options',   opts);

% 2a) MultiStart with more random starts
ms = MultiStart( ...
  'UseParallel',           true, ...
  'StartPointsToRun',      'bounds');

% nStarts = 40;  % try 40 different starts
% [estP,fval] = run(ms, problem, nStarts);

% ---- or ----

% 2b) Use GlobalSearch
gs = GlobalSearch();
[estP,fval] = run(gs, problem);


%% 5) Plot true vs recovered
paramNames = {'\mu_A','\mu_V','\lambda_A','\lambda_V','w'};
figure('Color','w','Position',[100 100 900 500]);
for i = 1:5
    subplot(2,3,i);
    scatter(trueP(i)*ones(nRecov,1), recovered(:,i), 36,'filled');
    hold on;
    plot(trueP(i), trueP(i),'r*','MarkerSize',12);
    xlabel('True'); ylabel('Recovered');
    title(paramNames{i},'Interpreter','tex');
    axis square; xl = xlim; ylim(xl);
    plot(xl,xl,'k:');
end

%% ———————————————————————————————————————————————
function SSD = computeSSD(params, lags, empMean, empSD, simN)
    % single‐call summary + SSD (no bootstrap)
    [pm, ps] = getSummaryRelay(params, lags, simN);
    SSD = sum( (pm - empMean).^2 + (ps - empSD).^2 );
end

function [predMean, predSD] = getSummaryRelay(params, lags, simN)
    % vectorized relay‐model: returns mean & SD at each lag
    muA = params(1); muV = params(2);
    lamA = params(3); lamV = params(4);
    w   = params(5);

    nL = numel(lags);
    predMean = zeros(1,nL);
    predSD   = zeros(1,nL);

    % loop is small (≤5), so this is fast
    for i=1:nL
        rt = getRelayLagRND(simN,muA,muV,lamA,lamV, ...
                            w,1-w,w,1-w,lags(i),false);
        rt = sampleDown(rt, 1000);
        predMean(i) = mean(rt);
        predSD(i)   = std(rt);
    end
end

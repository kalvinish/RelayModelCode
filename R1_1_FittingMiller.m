%% Miller (1982) Data Fitting with Inverse Gaussian Relay Model
% This script loads multisensory response time data, fits inverse Gaussian
% models to unisensory distributions, computes predicted multisensory
% CDFs under different models (Raab, Miller, relay), and compares them to
% empirical multisensory data.

clear; clc; close all;

% Add custom functions to MATLAB path
addpath(genpath(fullfile(pwd, 'Functions')));

%% Load Empirical Data
% Data files: '1.csv', '2.csv', '3.csv' correspond to A, V, and AV.
loadDir = fullfile(pwd, 'EmpiricalData', 'Miller82');
numQuantiles = 10;
empiricalData = nan(numQuantiles, 3);
probLevels = linspace(0.05, 0.95, numQuantiles);
loadedData = nan(numQuantiles, 2, 3);
fitParams = nan(3, 2);  % [mu, lambda] for each condition

for cond = 1:3
    % Read RT quantiles from CSV
    filePath = fullfile(loadDir, sprintf('%d.csv', cond));
    raw = readmatrix(filePath);
    empiricalData(:, cond) = raw(:, 1);

    % Build CDF data: [RT, probability]
    loadedData(:,1,cond) = raw(:,1);
    loadedData(:,2,cond) = probLevels;

    % Fit inverse Gaussian to unisensory quantile data
    [mu_opt, lambda_opt] = fitIG_fromCDF(loadedData(:,:,cond));
    fitParams(cond,:) = [mu_opt, lambda_opt];
end

%% Create Inverse Gaussian Distributions for A and V
xx = linspace(100, 800, 100);  % RT range for plotting
A_dist = makedist('InverseGaussian', 'mu', fitParams(1,1), 'lambda', fitParams(1,2));
V_dist = makedist('InverseGaussian', 'mu', fitParams(2,1), 'lambda', fitParams(2,2));

% Compute unisensory CDFs
cdf_A = cdf(A_dist, xx);
cdf_V = cdf(V_dist, xx);

%% Compute Race and RMI
cdf_Raab   = getRaabCDF(xx, A_dist.mu, V_dist.mu, A_dist.lambda, V_dist.lambda);
cdf_Miller = getMillerCDF(xx, A_dist.mu, V_dist.mu, A_dist.lambda, V_dist.lambda);

%% Optimal Relay Model Weight (RT-share) for AV Data
AV_data = loadedData(:,:,3);
[opt_w, fval] = getRTshare(AV_data, A_dist.mu, V_dist.mu, A_dist.lambda, V_dist.lambda);

% Relay CDF with optimal weight
cdf_AV = getRelayCDF(xx, A_dist.mu, V_dist.mu, A_dist.lambda, V_dist.lambda, opt_w, 1-opt_w, opt_w, 1-opt_w);
weights.optimal = opt_w;
weights.fval = fval;

%% Save Parameter Values
W.optimal_w = opt_w;
W.fval = fval;
params = {A_dist, V_dist, W};
save(fullfile(pwd, 'FittedParameters', 'Miller82_Parameters.mat'), 'params');

%% Compute Real-Stage Statistics for Relay Model
% First and second stage means and SDs for A and V channels
for ch = {'A', 'V'}
    dist = eval([ch{1}, '_dist']);
    w = opt_w;
    mu = dist.mu; lambda = dist.lambda;

    stats.(ch{1}).firstMean  = mu * w;
    stats.(ch{1}).firstSD    = sqrt((stats.(ch{1}).firstMean^3) / (lambda * w));
    stats.(ch{1}).secondMean = mu * (1 - w);
    stats.(ch{1}).secondSD   = sqrt((stats.(ch{1}).secondMean^3) / (lambda * (1 - w)));
end

%% Plot Fitted CDFs vs. Unisensory and Miller Predictions
figure; hold on;
plot(xx, cdf_A, 'g',  'LineWidth', 2);  % Auditory
plot(xx, cdf_V, 'b',  'LineWidth', 2);  % Visual
plot(xx, cdf_Raab, '--r', 'LineWidth', 2);  % Raab race
plot(xx, cdf_Miller, 'r', 'LineWidth', 2);  % Miller bound
plot(xx, cdf_AV, 'k',  'LineWidth', 2);  % Relay model

% Scatter empirical quantiles
cp = getCP(numQuantiles);
scatter(empiricalData(:,1), cp, 50, 'g',  'filled');
scatter(empiricalData(:,2), cp, 50, 'b',  'filled');
scatter(empiricalData(:,3), cp, 50, 'k',  'filled');

xlabel('Response Time (ms)'); ylabel('Cumulative Probability');
legend('A', 'V', 'Raab', 'Miller', 'Relay', 'Location', 'Southeast');
title('CDF Fits to Miller et al. (1982) Data');

%% Compute RSE and Violations
% Empirical gain and violations
rse_emp  = getGain(empiricalData);
viol_emp = getViolation(empiricalData);

% Model-based gain and violations from CDFs
griceCDF   = getGriceCDF(cdf_A, cdf_V);
model_rse  = getRSE_fromCDF(xx, cdf_AV, griceCDF);
model_viol = getViolation_fromCDF(xx, cdf_AV, cdf_Miller);

disp('Empirical vs. Model Comparison:');
disp(table(rse_emp, model_rse, viol_emp, model_viol, ...
    'VariableNames', {'RSE_Empirical','RSE_Model','Viol_Empirical','Viol_Model'}));

%% Evaluate Relay vs. Simple IG CDF Fit (R^2)
emp_IG = fitdist(empiricalData(:,3), 'InverseGaussian');
emp_IG_cdf = cdf(emp_IG, empiricalData(:,3));
predicted_CDF = getRelayCDF(empiricalData(:,3), A_dist.mu, V_dist.mu, A_dist.lambda, V_dist.lambda, opt_w, 1-opt_w, opt_w, 1-opt_w);

dataCDF = linspace(0.05, 0.95, numQuantiles);
R2_relay = corr(dataCDF', predicted_CDF).^2;
R2_ig    = corr(dataCDF', emp_IG_cdf).^2;

disp(['Relay model R^2 = ', num2str(R2_relay)]);
disp(['IG model    R^2 = ', num2str(R2_ig)]);

%% Final Combined CDF Plot with Shaded Areas (for Fig 1)
figure('Visible','on'); hold on;
xlim([100, 700]); yticks(0:0.2:1); box off;

% Shaded gain and violation regions
fillArea([empiricalData(:,3), getGrice(empiricalData(:,1:2))], [0.9 0.9 0.9], 1);
fillArea([empiricalData(:,3), getMiller(empiricalData(:,1:2))], [0.75 0.75 0.75], 1);

% Plot empirical points and model curves
markerOpts = {'o','MarkerFaceColor','w','MarkerSize',8,'LineStyle','-','LineWidth', 2};
plot(empiricalData(:,1), cp, markerOpts{:}, 'Color','g');
plot(empiricalData(:,2), cp, markerOpts{:}, 'Color','b');
plot(empiricalData(:,3), cp, markerOpts{:}, 'Color','k');
plot(getMiller(empiricalData(:,1:2)), cp, '-', 'Color','r', 'LineWidth', 2);
plot(getRaab(empiricalData(:,1:2)),  cp, '--','Color','r', 'LineWidth', 2);

xlabel('Response Time (ms)'); ylabel('Cumulative Probability');
legend({'','', 'A', 'V', 'AV', 'Miller Bound', 'Race Model'}, 'Location','Southeast','Box','off');

% Save figure
outFile = fullfile(pwd, 'Figures', 'miller82CDF.pdf');
exportgraphics(gcf, outFile, 'ContentType', 'vector');

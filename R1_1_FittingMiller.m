%% Fitting Miller 1982 Data with Inverse Gaussian

clear
clc
close all

addpath(genpath(fullfile('Functions')));

%% LOAD DATA TO WORKSPACE

loadPath = fullfile(cd, "EmpiricalData", "Miller82");

data = nan(10, 3);
params = nan(3,2);

for i = 1:3
    temp_data = readmatrix(fullfile(loadPath, [int2str(i), '.csv']));

    data(:,i) = temp_data(:,1);

    loadedData(:,1,i) = temp_data(:,1);
    loadedData(:,2,i) = linspace(0.05, 0.95, 10);

    % [mu_opt, lambda_opt, resnorm(i)] = fit_inverse_gaussian_cdf(loadedData(:,:,i));
    [mu_opt, lambda_opt, resnorm(i), ks_stat(i), h(i), p(i), boot_pvalue(i)] = fit_inverse_gaussian_cdf(loadedData(:,:,i));
    params(i,:) = [mu_opt, lambda_opt];
end

%% FIT INVERSE GAUSSIANS TO UNISENSORY

xx = linspace(100, 800, 100);

A = makedist("InverseGaussian", "mu", params(1,1), "lambda", params(1,2));
V = makedist("InverseGaussian", "mu", params(2,1), "lambda", params(2,2));

cdfA = cdf(A, xx);
cdfV = cdf(V, xx);

cdfRaab = myRaabCDF(xx, A.mu, V.mu, A.lambda, V.lambda);
cdfMiller = millerCDF(xx, A.mu, V.mu, A.lambda, V.lambda);

%% FIND OPTIMAL WEIGHT FOR MILLER 82 DATA

AVdata = readmatrix(fullfile(loadPath, [int2str(3), '.csv']));

AVdata(:,1) = AVdata(:,1);
AVdata(:,2) = linspace(0.05, 0.95, 10);

[W.optimal_w, W.fval, W.ks_stat, W.h, W.p, W.bootP] = find_optimal_weight(AVdata, A.mu, V.mu, A.lambda, V.lambda);

cdfAV = multiCDF(xx, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w);

params = {A, V, W};

realA.firstStageMean = A.mu * W.optimal_w;
realA.firstStageSD = sqrt((realA.firstStageMean^3)/(A.lambda * W.optimal_w));
realA.secondStageMean = A.mu * (1-W.optimal_w);
realA.secondStageSD = sqrt((realA.secondStageMean^3)/(A.lambda * (1-W.optimal_w)));

realV.firstStageMean = V.mu * W.optimal_w;
realV.firstStageSD = sqrt((realV.firstStageMean^3)/(V.lambda * W.optimal_w));
realV.secondStageMean = V.mu * (1-W.optimal_w);
realV.secondStageSD = sqrt((realV.secondStageMean^3)/(V.lambda * (1-W.optimal_w)));

% save(fullfile(cd, "FittedParameters", "Miller82_Parameters.mat"), "params")

%% PLOT FITTED ONTO EXTRACTED DATA

cp = getCP(size(data,1));

f1 = figure();

plot(xx, cdfA, 'g', LineWidth=2)
hold on
plot(xx, cdfV, 'b', LineWidth=2)
plot(xx, cdfRaab, '--r', LineWidth=2)
plot(xx, cdfMiller, 'r', LineWidth=2)

% Compute predicted CDFs at optimal weight and CI bounds
predicted_F_opt = multiCDF(xx, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w);

% Plot the predicted CDF at optimal weight
plot(xx, predicted_F_opt, 'k', 'LineWidth', 2);

scatter(data(:,1), cp, 'g', 'filled')
scatter(data(:,2), cp, 'b', 'filled')
scatter(data(:,3), cp, 'k', 'filled')

%% GET RSE AND VIOLATIONS

rse = getGain(data);
violation = getViolation(data);

miller = getMiller(data(:,1:2));
raab = getRaab(data(:,1:2));
grice = getGrice(data(:,1:2));

%% Get RSE and Violations from fitted parameters

griceCDF = getGriceCDF(uniCDF(xx, A.mu, A.lambda), uniCDF(xx, V.mu, V.lambda));
millercdf = millerCDF(xx, A.mu, V.mu, A.lambda, V.lambda);

% multi_rnd = sampleDown(multiRND(100000, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w),10);
% 
% scatter(multi_rnd, cp, 'r', 'filled')
% 
% my_rse = getGain([data(:,1:2) multi_rnd]);
% my_viol = getViolation([data(:,1:2) multi_rnd]);
% 
% disp([my_rse, rse])
% disp([my_viol, violation])

my_rse = getGainFromCDF(xx, predicted_F_opt, griceCDF);
my_viol = getViolationFromCDF(xx, predicted_F_opt, millercdf);

disp([my_rse, rse])
disp([my_viol, violation])

%% Predicted vs Empirical CDF

emp_data = data(:,3);

predicted_CDF = multiCDF(emp_data, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w);
% emp_norm = fitdist(emp_data, 'Normal');
emp_norm = fitdist(emp_data, 'InverseGaussian');
norm_cdf = cdf(emp_norm, emp_data);

empirical_CDF = linspace(0.05, 0.95, 10);

f2 = figure();
hold on;

plot(predicted_CDF, empirical_CDF, 'o');
plot([0,1],[0,1],'r--');   % 45° reference line
xlabel('Fitted CDF');
ylabel('Empirical CDF');
title('Empirical vs.\ Fitted CDF');

R2_relay = corr(empirical_CDF', predicted_CDF)^2;

f3 = figure();
hold on;

plot(norm_cdf, empirical_CDF, 'o');
plot([0,1],[0,1],'r--');   % 45° reference line
xlabel('Fitted CDF');
ylabel('Empirical CDF');
title('Empirical vs.\ Fitted CDF');

R2_ig = corr(empirical_CDF', norm_cdf)^2;

%% PLOT FINAL CDF

f4 = figure('Visible', 'on');
hold on
xmin = 100;
xmax = 700;

dataMrkrSize = 2.5;
mrkrSize = 8;
modelLineSize = 2.5;

fillArea([data(:,3) grice], [0.9 0.9 0.9], 1);
fillArea([data(:,3) miller], [0.75 0.75 0.75], 1);
plot(data(:,1), cp, "marker", "o", "MarkerSize", mrkrSize, "MarkerFaceColor", "w", "MarkerEdgeColor", "g", "Color", "g", "LineWidth", dataMrkrSize)
plot(data(:,2), cp, "marker", "o", "MarkerSize", mrkrSize, "MarkerFaceColor", "w", "MarkerEdgeColor", "b", "Color", "b", "LineWidth", dataMrkrSize)
plot(data(:,3), cp, "marker", "o", "MarkerSize", mrkrSize, "MarkerFaceColor", "w", "MarkerEdgeColor", "k", "Color", "k", "LineWidth", dataMrkrSize)
plot(miller, cp, "LineStyle", "-", "Color", "r", "LineWidth", modelLineSize)
plot(raab, cp, "LineStyle", "--", "Color", "r", "LineWidth", modelLineSize)

yticks(linspace(0, 1, 6))
xlim([xmin xmax])
xticks(linspace(xmin, xmax, 3))
box off
ax = gca;
ax.TickDir = "Out";
ax.YColor = "k";
ax.XColor = "k";
ax.LineWidth = 1.8;
ax.FontSize = 12;

ylabel("Cumulative Probability", "FontSize", 14, "Color", "k", "FontWeight", "bold")
xlabel("Response Time (s)", "FontSize", 14, "Color", "k", "FontWeight", "bold")

legendTxt = {"", "", "A", "V", "AV", "Miller's Bound", "Independent Race"};
legend(legendTxt, "location", "southeast", "box", "off")

% exportgraphics(f2, fullfile(cd, "Figures", "miller82CDF.pdf"), "ContentType", "vector")

%% SAMPLE 4 RANDOM TRIALS FOR FIGURE 1 TABLE

meanDiff = 0;

while (meanDiff < 65) || (meanDiff > 67)

    sample( :,1 ) = round( random( A, [4 1] ) );
    sample( :,2 ) = round( random( V, [4 1] ) );
    sample( :,3 ) = min( sample( :,1:2 ), [], 2 );
    
    descSample(1,:) = mean(sample);
    descSample(2,:) = std(sample);
    
    sampleRSE = getGain( sample );
    meanDiff = min( descSample( 1,1:2 ) ) - descSample( 1,3 );

end

%% Get LogL of relay and race

% [logL_relay, pdf_relay] = relayLogL(xx, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w);
% [logL_race, pdf_race] = raabLogL(xx, A.mu, V.mu, A.lambda, V.lambda);
% rt = multiRND(10000, A.mu, V.mu, A.lambda, V.lambda, W.optimal_w, 1-W.optimal_w, W.optimal_w, 1-W.optimal_w, 1);
% 
% f3 = figure();
% hold on
% 
% histogram(rt, "Normalization", "pdf")
% plot(xx, pdf_relay)
% plot(xx, pdf_race)

%% FUNCTIONS 

% function [optimal_w, fval, ks_stat, h, p] = find_optimal_weight(data, aMU, vMU, aLAMBDA, vLAMBDA)
%     % This function finds the optimal weight 'w' that minimizes the normalized RMSE
%     % between the empirical CDF of 'data' and the predicted CDF generated by
%     % 'multiCDF' using the given parameters.
%     %
%     % Inputs:
%     %   data     - N x 2 array where data(:,1) are the x values and data(:,2) are the empirical CDF values.
%     %   aMU, vMU, aLAMBDA, vLAMBDA - Parameters for the multiCDF function.
%     %
%     % Outputs:
%     %   optimal_w - Optimized weight 'w'.
%     %   fval      - Objective function value at the optimal weight.
%     %   ci_w      - 95% confidence interval for the optimized weight.
%     %   ks_stat   - Kolmogorov-Smirnov statistic for goodness-of-fit.
% 
%     empirical_F = data(:,2);
%     xx = data(:,1);
% 
%     % Remove duplicate xx values to ensure unique evaluation points
%     [xx, idx] = unique(xx);
%     empirical_F = empirical_F(idx);
% 
%     % Define the objective function to minimize
%     obj_fun = @(w) compute_distance(w, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA);
% 
%     % Define the bounds for the weight
%     lb = 0;
%     ub = 0.5;
% 
%     % Use fminbnd for optimization
%     options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
%     [optimal_w, fval] = fminbnd(obj_fun, lb, ub, options);
% 
%     % Compute the predicted CDF at optimal weight
%     predicted_F_opt = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, optimal_w, 1 - optimal_w, optimal_w, 1 - optimal_w);
% 
%     % Kolmogorov-Smirnov statistic
%     [h, p, ks_stat] = kstest(data(:,1), 'CDF', [xx, predicted_F_opt]);
% end

function [optimal_w, fval, ks_stat, h, p, boot_pvalue] = find_optimal_weight(data, aMU, vMU, aLAMBDA, vLAMBDA)
    % This function finds the optimal weight 'w' that minimizes the normalized RMSE
    % between the empirical CDF of 'data' and the predicted CDF generated by
    % 'multiCDF' using the given parameters. It then computes a KS statistic 
    % and uses bootstrapping to obtain an empirical p-value for the goodness-of-fit.
    %
    % Inputs:
    %   data      - N x 2 array where data(:,1) are the x values (raw observations)
    %               and data(:,2) are the empirical CDF values.
    %   aMU, vMU, aLAMBDA, vLAMBDA - Parameters for the multiCDF function.
    %
    % Outputs:
    %   optimal_w   - Optimized weight 'w'.
    %   fval        - Objective function value at the optimal weight.
    %   ks_stat     - KS statistic (maximum difference between the empirical CDF and the
    %                 fitted CDF).
    %   h           - KS test decision (1 if the null hypothesis is rejected at 5%
    %                 significance, 0 otherwise).
    %   boot_pvalue - Empirical p-value from the bootstrapping approach.
    
    % Extract the data: raw observations and empirical CDF
    empirical_F = data(:,2);
    xx = data(:,1);

    % Remove duplicate x-values (ensures unique, sorted evaluation points)
    [xx, idx] = unique(xx);
    empirical_F = empirical_F(idx);

    % Define the objective function to minimize (e.g., normalized RMSE between CDFs)
    obj_fun = @(w) compute_distance(w, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA);

    % Define the bounds for the weight
    lb = 0;
    ub = 0.5;

    % Optimize using fminbnd
    options = optimset('Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-10);
    [optimal_w, fval] = fminbnd(obj_fun, lb, ub, options);

    % Compute the predicted (fitted) CDF at the optimal weight using multiCDF.
    % The arguments for multiCDF below are as per your previous call.
    predicted_F_opt = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, ...
                               optimal_w, 1 - optimal_w, optimal_w, 1 - optimal_w);

    % Compute the KS statistic using the raw data and the fitted CDF.
    % Ensure that the hypothesized CDF matrix spans the data range.
    [h, p, ks_stat] = kstest(data(:,1), 'CDF', [xx, predicted_F_opt]);

    % --- Bootstrapping to approximate the null distribution ---
    % Number of bootstrap iterations
    B = 10000;
    n = length(data(:,1));  % sample size
    ks_bootstrap = zeros(B,1);

    % Create an inverse CDF function using interpolation.
    % This assumes that predicted_F_opt is monotonic and spans from 0 to 1.
    inv_cdf = @(u) interp1(predicted_F_opt, xx, u, 'linear', 'extrap');

    for b = 1:B
        % Generate a bootstrap sample:
        % Draw n independent uniform random numbers
        u = rand(n, 1);
        % Generate synthetic sample using inverse transform sampling.
        x_sim = inv_cdf(u);
        
        % Ensure that x_sim values lie within the range of xx to satisfy kstest requirements.
        x_sim = min(max(x_sim, min(xx)), max(xx));
        
        % Compute the KS statistic for the simulated sample versus the fitted CDF.
        % The hypothesized CDF [xx, predicted_F_opt] now spans the sample.
        [~, ~, ks_bootstrap(b)] = kstest(x_sim, 'CDF', [xx, predicted_F_opt]);
    end

    % Compute the empirical (bootstrap) p-value:
    % p-value = proportion of bootstrap KS statistics that exceed the observed ks_stat.
    boot_pvalue = mean(ks_bootstrap >= ks_stat);
end




%% Subfunction to compute the RMSE between CDFs

function dist = compute_distance(w, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA)
    % Compute predicted CDF using multiCDF
    predicted_F = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, w, 1 - w, w, 1 - w);

    % Compute residuals
    residuals = empirical_F - predicted_F;
    
    % Compute root mean square error (RMSE)
    dist = sqrt(mean(residuals.^2));
end

%%

% function [mu_opt, lambda_opt, resnorm, ks_stat, h, p, boot_pvalue] = fit_inverse_gaussian_cdf(data)
%     % This function fits an Inverse Gaussian distribution to an empirical CDF.
%     % It returns the estimated parameters (mu and lambda), the residual norm,
%     % the KS statistic comparing the fitted CDF to the data, and an empirical 
%     % p-value (via bootstrapping) for the goodness-of-fit test.
%     %
%     % Input:
%     %   data - an N x 2 array where data(:,1) are the x values and data(:,2) are
%     %          the empirical CDF values at those x values.
%     %
%     % Outputs:
%     %   mu_opt      - Estimated mu (mean) parameter of the Inverse Gaussian distribution
%     %   lambda_opt  - Estimated lambda (shape) parameter of the Inverse Gaussian distribution
%     %   resnorm     - Residual norm from the least-squares fit
%     %   ks_stat     - KS statistic (maximum absolute difference between the empirical CDF and
%     %                 the fitted CDF)
%     %   h           - KS test decision (1 if the null hypothesis is rejected at 5% significance level,
%     %                 0 otherwise)
%     %   p           - p-value from the KS test based on the standard asymptotic distribution
%     %   boot_pvalue - Empirical p-value obtained from bootstrapping
%     %
%     % Rationale:
%     % The KS test is used here as a measure of the maximum discrepancy between the observed and fitted CDFs
%     % [Massey, 1951]. However, since parameters are estimated from the data, the asymptotic null distribution is 
%     % not strictly valid [Lilliefors, 1967]. To address this, we use bootstrapping [Efron, 1979] by generating 
%     % synthetic samples using the inverse CDF (quantile function) and comparing their KS statistics to that observed.
% 
%     % Extract x values and empirical CDF values from the data
%     x = data(:,1);
%     F_empirical = data(:,2);
% 
%     % Remove duplicate x values for a unique, sorted evaluation grid
%     [x, idx] = unique(x);
%     F_empirical = F_empirical(idx);
% 
%     % Initial guesses for the parameters
%     mu0 = mean(x);
%     lambda0 = 1;
%     params0 = [mu0 lambda0];
% 
%     % Define the model function using MATLAB's built-in Inverse Gaussian CDF
%     model_fun = @(params, x) cdf('InverseGaussian', x, params(1), params(2));
% 
%     % Define the residual function to be minimized (difference between empirical and model CDF)
%     function residuals = simple_residuals(params)
%         F_model = model_fun(params, x);
%         residuals = F_empirical - F_model;
%     end
% 
%     % Optimization options for lsqnonlin
%     options = optimoptions('lsqnonlin', 'Display', 'off', ...
%                            'TolFun', 1e-10, 'TolX', 1e-10);
% 
%     % Run the least-squares optimization to fit the Inverse Gaussian CDF
%     [params_opt, resnorm] = lsqnonlin(@simple_residuals, params0, [], [], options);
% 
%     % Extract optimized parameters
%     mu_opt = params_opt(1);
%     lambda_opt = params_opt(2);
% 
%     % Compute the fitted CDF at the evaluation points x
%     F_fit = cdf('InverseGaussian', x, mu_opt, lambda_opt);
% 
%     % Perform the standard KS test using the raw observations (x)
%     % The hypothesized CDF is provided as a two-column matrix [x, F_fit].
%     [h, p, ks_stat] = kstest(x, 'CDF', [x, F_fit]);
% 
%     % --- Bootstrapping to adjust for parameter estimation uncertainty ---
%     % Set the number of bootstrap iterations
%     B = 10000;
%     n = length(x);
%     ks_bootstrap = zeros(B,1);
% 
%     % Generate bootstrap samples using the inverse CDF (quantile function)
%     % The inverse CDF is computed using MATLAB's icdf function for the Inverse Gaussian distribution.
%     for b = 1:B
%         u = rand(n, 1);  % generate n uniform random numbers in [0,1]
%         x_sim = icdf('InverseGaussian', u, mu_opt, lambda_opt);
%         % Ensure simulated values fall within the interval spanned by x
%         x_sim = min(max(x_sim, min(x)), max(x));
%         [~, ~, ks_bootstrap(b)] = kstest(x_sim, 'CDF', [x, F_fit]);
%     end
% 
%     % The empirical (bootstrap) p-value is computed as the proportion of bootstrap
%     % KS statistics exceeding the observed ks_stat.
%     boot_pvalue = mean(ks_bootstrap >= ks_stat);
% end

function [mu_opt, lambda_opt, rmse, ks_stat, h, p, boot_pvalue] = fit_inverse_gaussian_cdf(data)
    % Fits an Inverse Gaussian CDF by directly minimising RMSE.
    %
    % Inputs:
    %   data - N×2: [x, emp_CDF(x)]
    %
    % Outputs:
    %   mu_opt, lambda_opt — fitted parameters
    %   rmse      — root-mean-square error of the fit
    %   ks_stat,h,p — KS test results vs. fitted CDF
    %   boot_pvalue — bootstrap KS p-value

    % Extract and dedupe
    x           = data(:,1);
    F_empirical = data(:,2);
    [x, idx]    = unique(x);
    F_empirical = F_empirical(idx);
    n           = numel(x);

    % Initial guess
    mu0     = mean(x);
    lambda0 = 1;
    params0 = [mu0, lambda0];

    % Model CDF handle
    modelCDF = @(params) cdf('InverseGaussian', x, params(1), params(2));

    % Objective: RMSE
    objfun = @(params) sqrt( mean( (F_empirical - modelCDF(params)).^2 ) );

    % Constrain mu>0, lambda>0
    lb = [eps, eps];
    ub = [];

    opts = optimoptions('fmincon', ...
               'Display','off', ...
               'TolFun',1e-10, ...
               'TolX',1e-10);

    % Run the minimisation
    [params_opt, rmse] = fmincon(objfun, params0, [],[],[],[], lb, ub, [], opts);

    mu_opt     = params_opt(1);
    lambda_opt = params_opt(2);

    % Compute fitted CDF and do KS test
    F_fit       = modelCDF(params_opt);
    [h, p, ks_stat] = kstest(x, 'CDF', [x, F_fit]);

    % Bootstrap for KS p-value
    B             = 10000;
    ks_bootstrap  = zeros(B,1);
    for b = 1:B
        u       = rand(n,1);
        x_sim   = icdf('InverseGaussian', u, mu_opt, lambda_opt);
        x_sim   = min(max(x_sim, min(x)), max(x));
        [~, ~, ks_bootstrap(b)] = kstest(x_sim, 'CDF', [x, F_fit]);
    end
    boot_pvalue = mean(ks_bootstrap >= ks_stat);
end



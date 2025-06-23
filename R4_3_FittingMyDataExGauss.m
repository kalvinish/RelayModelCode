%% Fitting Miller 1982 Data with ExGaussian

clear
clc
close all

addpath(genpath(fullfile('Functions')));

%% Load in empirical data from two button exp

loadPath = fullfile(cd, "EmpiricalData", "sampledRTs_4Relay.mat");
data = load(loadPath).save_samples;

n = size(data,1);
cp = getCP(n);

for ptcpt = 1:size(data,3)
    for cond = 1:2
        [mu_opt, sigma_opt, tau_opt] = fit_exgaussian_cdf([data(:,cond,ptcpt),cp]);
        params(cond,:,ptcpt) = [mu_opt, sigma_opt, tau_opt];
    end
end

% Plot fits
n = 1000;
xMin = 0.1;
xMax = 0.7;
xx = linspace(xMin, xMax, n);

f1 = figure();
hold on

vincent_data = mean(data, 3);

scatter(vincent_data(:,1), cp, 30, 'g', 'filled')
scatter(vincent_data(:,2), cp, 30, 'b', 'filled')
scatter(vincent_data(:,3), cp, 30, 'r', 'filled')

for ptcpt = 1:size(data,3)
    rr(:,1,ptcpt) = sort(exgauss_rnd(params(1,1,ptcpt), params(1,2,ptcpt), params(1,3,ptcpt), n, 1));
    rr(:,2,ptcpt) = sort(exgauss_rnd(params(2,1,ptcpt), params(2,2,ptcpt), params(2,3,ptcpt), n, 1));
    rr(:,3,ptcpt) = exgauss_relay_rnd(params(:,1,ptcpt), params(:,2,ptcpt), params(:,3,ptcpt), n, 1);
    rr(:,4,ptcpt) = sort(exgauss_race_rnd(params(:,1,ptcpt), params(:,2,ptcpt), params(:,3,ptcpt), n));
end

vincent_rr = mean(rr,3);
cp = getCP(n);
 
plot(vincent_rr(:,1), cp, 'Color', 'g', 'LineWidth', 2)
plot(vincent_rr(:,2), cp, 'Color', 'b', 'LineWidth', 2)
plot(vincent_rr(:,3), cp, 'Color', 'r', 'LineWidth', 2)
plot(vincent_rr(:,4), cp, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--')

xlim([xMin xMax])

%% Functions

function [mu_opt, sigma_opt, tau_opt] = fit_exgaussian_cdf(data)
    % This function fits an ExGaussian distribution to an empirical CDF.
    % It returns the estimated parameters mu, sigma and tau.
    %
    % Input:
    %   data - an N x 2 array where data(:,1) are the x values and data(:,2) are
    %          the empirical CDF values at those x values.
    %
    % Output:
    %   mu_opt     - Estimated mu parameter of the Gaussian component
    %   sigma_opt  - Estimated sigma parameter of the Gaussian component
    %   tau_opt    - Estimated tau parameter of the Exponential component

    x = data(:,1);
    F_empirical = data(:,2);

    % Remove duplicate x values to ensure unique evaluation points
    [x, idx] = unique(x);
    F_empirical = F_empirical(idx);

    % Initial guesses for the parameters
    mu0 = mean(x);
    sigma0 = 0.1;
    tau0 = 0.1;

    params0 = [mu0 sigma0 tau0];

    % Define the model function using built-in Inverse Gaussian CDF
    model_fun = @(params, x) exgauss_cdf(x, params(1), params(2), params(3));

    % Define the weighted residual function
    function residuals = simple_residuals(params)
        F_model = model_fun(params, x);
        residuals = F_empirical - F_model;
    end

    % Optimization options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'TolFun', 1e-10, 'TolX', 1e-10);

    % Run the optimization
    [params_opt, ~] = lsqnonlin(@simple_residuals, params0, [], [], options);

    % Extract the optimized parameters
    mu_opt = params_opt(1);
    sigma_opt = params_opt(2);
    tau_opt = params_opt(3);
end

% Race between two exponential distributions (RND)
function rr = exponential_race_rnd(tau, n)
    rr1 = exprnd(tau(1), n, 1);
    rr2 = exprnd(tau(2), n, 1);

    rr = min([rr1 rr2], [], 2);
end

% Race between two Gaussian distributions (RND)
function rr = gaussian_race_rnd(mu, sigma, n)
    rr1 = normrnd(mu(1), sigma(1), n, 1);
    rr2 = normrnd(mu(2), sigma(2), n, 1);

    rr = min([rr1 rr2], [], 2);
end

function rr = exgauss_race_rnd(mu, sigma, tau, n)
    rr1 = exgauss_rnd(mu(1), sigma(1), tau(1), n, 1);
    rr2 = exgauss_rnd(mu(2), sigma(2), tau(2), n, 1);

    rr = min([rr1 rr2], [], 2);
end

function rr = exgauss_relay_rnd(mu, sigma, tau, n, sort_flg)
    rr1 = exponential_race_rnd(tau, n);
    rr2 = gaussian_race_rnd(mu, sigma, n);

    if ~exist('sort_flg','var')
        rr = rr1 + rr2;
    elseif sort_flg == 1
        rr = sort(rr1 + rr2);
    else
        rr = rr1 + rr2;
    end
end
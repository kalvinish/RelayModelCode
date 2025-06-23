%% Fitting Miller 1982 Data with ExGaussian

clear
clc
close all

addpath(genpath(fullfile('Functions')));

%% Load data and fit unisensory conditions with exgaussian

loadPath = fullfile(cd, "EmpiricalData", "Miller82");

data = nan(10, 3);
params = nan(3,3);

for i = 1:3
    temp_data = readmatrix(fullfile(loadPath, [int2str(i), '.csv']));

    data(:,i) = temp_data(:,1) / 1000;
    F_emp(:,i) = temp_data(:,2);

    loadedData(:,1,i) = temp_data(:,1) / 1000;
    loadedData(:,2,i) = linspace(0.05, 0.95, 10);

    [mu_opt, sigma_opt, tau_opt] = fit_exgaussian_cdf(loadedData(:,:,i));
    params(i,:) = [mu_opt, sigma_opt, tau_opt];
end

%% Plot fits

n = 1000;
xMin = 0.2;
xMax = 0.8;
xx = linspace(xMin, xMax, n);

f1 = figure();
hold on
my_linewidth = 2;
my_markersize = 60;

cp = getCP(size(data,1));

scatter(data(:,1),cp,my_markersize,'g','filled','o','DisplayName','Emp-A')
scatter(data(:,2),cp,my_markersize,'b','filled','o','DisplayName','Emp-V')
scatter(data(:,3),cp,my_markersize,'r','filled','o','DisplayName','Emp-AV')

yy1 = exgauss_cdf(xx, params(1,1), params(1,2), params(1,3));
yy2 = exgauss_cdf(xx, params(2,1), params(2,2), params(2,3));
yy3 = exgauss_relay_cdf(xx, params(:,1), params(:,2), params(:,3));
race_yy = exgauss_race_cdf(xx, params(:,1), params(:,2), params(:,3));
miller_yy = exgauss_miller_cdf(xx, params(:,1), params(:,2), params(:,3));

plot(xx,yy1,'Color','g','LineWidth',my_linewidth,'DisplayName','Model-A')
plot(xx,yy2,'Color','b','LineWidth',my_linewidth,'DisplayName','Model-V')
plot(xx,yy3,'Color','r','LineWidth',my_linewidth,'DisplayName','Model-AV')
plot(xx,race_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','Race', 'LineStyle','--')
plot(xx,miller_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','RMI', 'LineStyle','-.')

legend('Location','southeast', 'Box','off')
xlim([xMin xMax])
ylim([0 1])
title('Model fitted to unisensory data only')

%% Fit to all data, including AV

params_opt = fit_exgaussian_cdf_full(data, F_emp);

% Plot fits
n = 1000;
xMin = 0.2;
xMax = 0.8;
xx = linspace(xMin, xMax, n);

f2 = figure();
hold on
my_linewidth = 2;
my_markersize = 60;

cp = getCP(size(data,1));

scatter(data(:,1),cp,my_markersize,'g','filled','o','DisplayName','Emp-A')
scatter(data(:,2),cp,my_markersize,'b','filled','o','DisplayName','Emp-V')
scatter(data(:,3),cp,my_markersize,'r','filled','o','DisplayName','Emp-AV')

yy1 = exgauss_cdf(xx, params_opt(1,1), params_opt(1,2), params_opt(1,3));
yy2 = exgauss_cdf(xx, params_opt(2,1), params_opt(2,2), params_opt(2,3));
yy3 = exgauss_relay_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));
race_yy = exgauss_race_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));
miller_yy = exgauss_miller_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));

plot(xx,yy1,'Color','g','LineWidth',my_linewidth,'DisplayName','Model-A')
plot(xx,yy2,'Color','b','LineWidth',my_linewidth,'DisplayName','Model-V')
plot(xx,yy3,'Color','r','LineWidth',my_linewidth,'DisplayName','Model-AV')
plot(xx,race_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','Race', 'LineStyle','--')
plot(xx,miller_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','RMI', 'LineStyle','-.')

legend('Location','southeast', 'Box','off')
xlim([xMin xMax])
ylim([0 1])
title('Model fitted to uni and multisensory data')

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

function params_opt = fit_exgaussian_cdf_full(data, F_emp)
    % Define the model function using built-in Inverse Gaussian CDF
    model_fun = @(params, x) exgauss_relay_full_cdf(x, params(:,1), params(:,2), params(:,3));
    
    params0 = [0.3 0.1 0.1;...
               0.3 0.1 0.1];
    
    % Define the weighted residual function
    function residuals = simple_residuals(params)
        F_model = model_fun(params, data);
        residuals = F_emp - F_model;
    end
    
    % Optimization options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'TolFun', 1e-10, 'TolX', 1e-10);
    
    % Run the optimization
    [params_opt, ~] = lsqnonlin(@simple_residuals, params0, [], [], options);
end

function yy = exgauss_relay_full_cdf(xx, mu, sigma, tau)
    yy1 = exgauss_cdf(xx(:,1), mu(1), sigma(1), tau(1));
    yy2 = exgauss_cdf(xx(:,2), mu(2), sigma(2), tau(2));
    yy3 = exgauss_relay_cdf(xx(:,3), mu, sigma, tau)';

    yy = [yy1 yy2 yy3];
end
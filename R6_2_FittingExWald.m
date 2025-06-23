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

    [mu_opt, sigma_opt, tau_opt] = fit_exwald_cdf(loadedData(:,:,i));
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

yy1 = exwald_cdf(xx, params(1,1), params(1,2), params(1,3));
yy2 = exwald_cdf(xx, params(2,1), params(2,2), params(2,3));
yy3 = exwald_relay_cdf(xx, params(:,1), params(:,2), params(:,3));
% race_yy = exgauss_race_cdf(xx, params(:,1), params(:,2), params(:,3));
miller_yy = exwald_miller_cdf(xx, params(:,1), params(:,2), params(:,3));

plot(xx,yy1,'Color','g','LineWidth',my_linewidth,'DisplayName','Model-A')
plot(xx,yy2,'Color','b','LineWidth',my_linewidth,'DisplayName','Model-V')
plot(xx,yy3,'Color','r','LineWidth',my_linewidth,'DisplayName','Model-AV')
% plot(xx,race_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','Race', 'LineStyle','--')
plot(xx,miller_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','RMI', 'LineStyle','-.')

legend('Location','southeast', 'Box','off')
xlim([xMin xMax])
ylim([0 1])
title('Model fitted to unisensory data only')

%% Fit to all data, including AV

params_opt = fit_exwald_cdf_full(data, F_emp);

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

yy1 = exwald_cdf(xx, params_opt(1,1), params_opt(1,2), params_opt(1,3));
yy2 = exwald_cdf(xx, params_opt(2,1), params_opt(2,2), params_opt(2,3));
yy3 = exwald_relay_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));
% race_yy = exgauss_race_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));
miller_yy = exwald_miller_cdf(xx, params_opt(:,1), params_opt(:,2), params_opt(:,3));

plot(xx,yy1,'Color','g','LineWidth',my_linewidth,'DisplayName','Model-A')
plot(xx,yy2,'Color','b','LineWidth',my_linewidth,'DisplayName','Model-V')
plot(xx,yy3,'Color','r','LineWidth',my_linewidth,'DisplayName','Model-AV')
% plot(xx,race_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','Race', 'LineStyle','--')
plot(xx,miller_yy,'Color','r','LineWidth',my_linewidth,'DisplayName','RMI', 'LineStyle','-.')

legend('Location','southeast', 'Box','off')
xlim([xMin xMax])
ylim([0 1])
title('Model fitted to uni and multisensory data')

%% Functions

function [mu_opt, lambda_opt, tau_opt] = fit_exwald_cdf(data)
    % Extract x & empirical F
    x = data(:,1);
    F_empirical = data(:,2);
    [x, idx] = unique(x);
    F_empirical = F_empirical(idx);

    % Initial guesses
    mu0     = mean(x);
    lambda0 = 4;
    tau0    = 0.1;
    params0 = [mu0, lambda0, tau0];

    % model and residual
    model_fun = @(p,x) exwald_cdf(x, p(1), p(2), p(3));
    resid     = @(p) F_empirical - model_fun(p,x);

    % bounds: only lambda (p(2)) must be >0
    % we use eps instead of 0 to prevent the solver pinning exactly at zero
    lb = [eps, eps,       -Inf];
    ub = [ Inf,  Inf,       Inf];

    opts = optimoptions('lsqnonlin','Display','off',...
                        'TolFun',1e-10,'TolX',1e-10);

    params_opt = lsqnonlin(resid, params0, lb, ub, opts);

    mu_opt     = params_opt(1);
    lambda_opt = params_opt(2);
    tau_opt    = params_opt(3);
end


function cdf = exwald_cdf(xx, mu, lambda, t, N)
%EXWALD_CDF_FAST  Fast CDF of the ex-Wald via a single cumtrapz
%
%   CDF = EXWALD_CDF_FAST(XX, MU, LAMBDA, T) does the same as
%   EXWALD_CDF, but builds one grid of N points over [0,max(XX)],
%   computes
%       integrand(u) = f_W(u)*exp(u/t)
%   once, and then uses cumtrapz+interp1 to get
%       I(x) = ∫₀ˣ f_W(u) e^{u/t} du
%   for all x in XX in one vectorized swoop.
%
%   You can optionally pass N (default 2 000) to control grid density.

    if nargin<5, N = 2000; end

    xx = xx(:);
    xmax = max(xx);
    if xmax<=0
        cdf = zeros(size(xx));
        return
    end

    % 1) Build the inverse‐Gaussian object
    pd = makedist('InverseGaussian','mu',mu,'lambda',lambda);

    % 2) Create a common grid u and evaluate integrand = f_W(u)*exp(u/t)
    u = linspace(0, xmax, N).';
    fw = pd.pdf(u);                   % inverse‐Gaussian density
    integ = fw .* exp(u./t);          % f_W(u)*e^{u/t}

    % 3) Do one cumulative‐integral
    cumI = cumtrapz(u, integ);       % cumI(k) ≈ ∫₀^{u(k)} f_W(u)e^{u/t} du

    % 4) For each xx(i), interpolate that partial‐integral
    Ixx = interp1(u, cumI, xx, 'linear', 0);

    % 5) Assemble the CDF vectorized
    cdf = pd.cdf(xx) - exp(-xx./t).*Ixx;
end

% Race between two exponential distributions (CDF)
function yy = exponential_race_cdf(xx, tau)
    F1 = expcdf(xx, tau(1));
    F2 = expcdf(xx, tau(2));

    yy = F1 + F2 - (F1 .* F2);
end

% Race between two wald distributions (PDF)
function yy = wald_race_pdf(xx, mu, lambda)

    pd1 = makedist('InverseGaussian','mu',mu(1),'lambda',lambda(1));
    pd2 = makedist('InverseGaussian','mu',mu(2),'lambda',lambda(2));

    f1 = pd1.pdf(xx);
    F1 = pd1.cdf(xx);
    f2 = pd2.pdf(xx);
    F2 = pd2.cdf(xx);
    
    yy = f1 .* (1 - F2) + f2 .* (1 - F1);
    
    yy(isnan(yy)) = 0;
end

% Relay between exgaussian distributions (CDF)
function yy = exwald_relay_cdf(xx, mu, lambda, tau)
    parfor ii=1:length(xx)
        if xx(ii)>0
            fun = @(t) exponential_race_cdf(xx(ii)-t, tau) .*...
                       wald_race_pdf(t, mu, lambda);

            yy(ii) = integral( fun, -Inf, Inf );
        end
    end 
end

% Race Model Inequality for ExGaussian
function yy = exwald_miller_cdf(xx, mu, lambda, tau)
    F1 = exwald_cdf(xx, mu(1), lambda(1), tau(1));
    F2 = exwald_cdf(xx, mu(2), lambda(2), tau(2));
    
    yy = F1 + F2;
end

function params_opt = fit_exwald_cdf_full(data, F_emp)
    % Define the model function using built-in Inverse Gaussian CDF
    model_fun = @(params, x) exwald_relay_full_cdf(x, params(:,1), params(:,2), params(:,3));
    
    params0 = [0.3 4 0.1;...
               0.3 4 0.1];
    
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

function yy = exwald_relay_full_cdf(xx, mu, lambda, tau)
    yy1 = exwald_cdf(xx(:,1), mu(1), lambda(1), tau(1));
    yy2 = exwald_cdf(xx(:,2), mu(2), lambda(2), tau(2));
    yy3 = exwald_relay_cdf(xx(:,3), mu, lambda, tau)';

    yy = [yy1 yy2 yy3];
end
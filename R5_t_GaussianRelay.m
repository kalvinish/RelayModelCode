%% GETTING EXGAUSSIAN RSE + VIOLATIONS
clear; clc; close all
addpath(genpath(fullfile('Functions')));

%% Parameters

% Weight parameter (fraction of the total mean/variance in stage 1)
w = 1;  

% Overall RT distribution parameters for auditory (A) and visual (V)
mu_total    = [0.4, 0.4];    % [A_overall_mu, V_overall_mu]
sigma_total = [0.05, 0.05];  % [A_overall_sigma, V_overall_sigma]

% Split into two serial stages per modality via helper function
[mu, sigma] = split_params(mu_total, sigma_total, w);

% Simulation settings
n    = 100000;
simN = 10000;
xMin = 0;
xMax = 1.1;
xx   = linspace(xMin, xMax, n);

%% Test unisensory

yy   = gauss_uni_cdf(xx, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));
rndA = gauss_uni_rnd(n,    mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));

f1 = figure();
tiledlayout(1,3,'TileSpacing','compact')

nexttile; hold on
plot(xx, normcdf(xx, mu_total(1), sigma_total(1)),'LineWidth',2,'DisplayName','Actual CDF')
plot(xx,yy,'.','LineWidth',2,'DisplayName','CDF');
plot(sort(rndA), getCP(n),'--','LineWidth',2,'DisplayName','Rnd');
legend('Box','off')
title('Gaussian Unisensory (A)')
xlabel('RT (s)'); ylabel('CDF')

%% Test gaussian race

yyRv = gauss_race_cdf(xx, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2));
rndRv = gauss_race_rnd(n, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2));

nexttile; hold on
plot(xx,yyRv,'LineWidth',2,'DisplayName','Race CDF');
plot(sort(rndRv),getCP(n),'--','LineWidth',2,'DisplayName','Race Rnd');
legend('Box','off')
title('Gaussian Race (A vs V)')
xlabel('RT (s)'); ylabel('CDF')

%% Test gaussian relay

yyRelay = gauss_relay_cdf(xx, mu, sigma);
rndRelay = gauss_relay_rnd(n, mu, sigma);

nexttile; hold on
plot(xx,yyRelay,'LineWidth',2,'DisplayName','Relay CDF');
plot(sort(rndRelay),getCP(n),'--','LineWidth',2,'DisplayName','Relay Rnd');
legend('Box','off')
title('Gaussian Relay')
xlabel('RT (s)'); ylabel('CDF')

%% Plot CDF of uni + multi conditions compared to RMI

aYY     = gauss_uni_cdf(xx, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));
vYY     = gauss_uni_cdf(xx, mu(1,2), sigma(1,2), mu(2,2), sigma(2,2));
relayYY = gauss_relay_cdf(xx, mu, sigma);
simXX   = gauss_relay_rnd(simN, mu, sigma);
simYY   = getCP(simN);
millerYY= gauss_miller_cdf(xx, mu, sigma);
griceYY = getGriceCDF(aYY, vYY);
violation = getViolationFromCDF(xx, relayYY, millerYY);
gain      = getGainFromCDF(xx, relayYY, griceYY);

f2 = figure();
tiledlayout(1,3)
lw = 2;

nexttile; hold on
plot(xx, aYY,     'LineWidth',lw,'DisplayName','A-only')
plot(xx, vYY,     'LineWidth',lw,'DisplayName','V-only')
plot(sort(simXX), simYY, 'LineWidth',lw,'DisplayName','sim AV')
plot(xx, relayYY, 'LineWidth',lw,'DisplayName','Relay AV')
plot(xx, millerYY,'LineWidth',lw,'DisplayName','RMI')
xlim([xMin, xMax]); ylim([0 1])
legend('Box','off','Location','southeast')
xlabel('RT (s)'); ylabel('CDF')

%% Get RSE and Violations across different weight values

wVals     = linspace(0,0.5,100);
violations = zeros(size(wVals));
gains      = zeros(size(wVals));

parfor i = 1:length(wVals)
    ww = wVals(i);
    [muL, sigmaL] = split_params(mu_total, sigma_total, ww);
    
    aYY_i      = gauss_uni_cdf(xx, muL(1,1), sigmaL(1,1), muL(2,1), sigmaL(2,1));
    vYY_i      = gauss_uni_cdf(xx, muL(1,2), sigmaL(1,2), muL(2,2), sigmaL(2,2));
    relayYY_i  = gauss_relay_cdf(xx, muL, sigmaL);
    millerYY_i = gauss_miller_cdf(xx, muL, sigmaL);
    griceYY_i  = getGriceCDF(aYY_i, vYY_i);

    violations(i) = getViolationFromCDF(xx, relayYY_i, millerYY_i);
    gains(i)      = getGainFromCDF(xx, relayYY_i, griceYY_i);
end

nexttile;
plot(wVals, gains,'-o')
xlabel('Weight w'); ylabel('RSE Gain');
title('RSE vs. split weight')

nexttile;
plot(wVals, violations,'-o')
xlabel('Weight w'); ylabel('Violation');
title('Violation vs. split weight')


%% FUNCTIONS
function [mu, sigma] = split_params(mu_total, sigma_total, w)
    % mu:     2×2 matrix of stage means
    % sigma:  2×2 matrix of stage SDs
    % that swap under w <-> 1-w

    % means still split linearly
    mu1 = w .* mu_total;
    mu2 = (1 - w) .* mu_total;

    % variances split linearly, so SDs split by sqrt(w)
    s1 = sqrt(w)   .* sigma_total;
    s2 = sqrt(1-w) .* sigma_total;

    mu    = [mu1;  mu2];
    sigma = [s1;   s2];
end

function yy = gauss_uni_cdf(x, mu1, sigma1, mu2, sigma2)
    mu_Z = mu1 + mu2;
    sigma_Z = sqrt(sigma1^2 + sigma2^2);

    % Compute the CDF using the error function
    yy = 0.5 * (1 + erf((x - mu_Z) / (sqrt(2) * sigma_Z)));
end

function rnd = gauss_uni_rnd(n, mu1, sigma1, mu2, sigma2)
    rnd1 = mu1 + randn(n,1) .* sigma1;
    rnd2 = mu2 + randn(n,1) .* sigma2;
    
    rnd = rnd1 + rnd2;
end

function yy = gauss_race_cdf(xx, mu1, sigma1, mu2, sigma2)
    F1 = normcdf(xx, mu1, sigma1);
    F2 = normcdf(xx, mu2, sigma2);

    yy = F1 + F2 - (F1 .* F2);
end

function rnd = gauss_race_rnd(n, mu1, sigma1, mu2, sigma2)
    rnd1 = mu1 + randn(n,1) .* sigma1;
    rnd2 = mu2 + randn(n,1) .* sigma2;

    rnd = min([rnd1 rnd2], [], 2);
end

% Race between two Gaussian distributions (PDF)
function yy = gaussian_race_pdf(xx, mu1, sigma1, mu2, sigma2)
    F1 = normcdf(xx, mu1, sigma1);
    F2 = normcdf(xx, mu2, sigma2);
    f1 = normpdf(xx, mu1, sigma1);
    f2 = normpdf(xx, mu2, sigma2);
    
    yy = f1 .* (1 - F2) + f2 .* (1 - F1);
    
    yy(isnan(yy)) = 0;
end

function yy = gauss_relay_cdf(xx, mu, sigma)
%GAUSS_RELAY_CDF serial relay CDF with zero‐sigma fallback
%   xx:    time vector
%   mu:    2×2 matrix of stage means [stage×modality]
%   sigma: 2×2 matrix of stage SDs    [stage×modality]

    %--- detect degenerate stage(s) ---
    if any(sigma(:) == 0)
        % rebuild overall mu & sigma for each modality
        mu_total    = sum(mu,    1);             % 1×2: [A_overall, V_overall]
        sigma_total = sqrt( sum(sigma.^2, 1) );  % 1×2: combine variances
        
        % fallback to simple race CDF across modalities
        yy = gauss_race_cdf(xx, ...
             mu_total(1), sigma_total(1), ...
             mu_total(2), sigma_total(2));
        return;
    end

    %--- normal relay case ---
    yy = zeros(size(xx));
    parfor ii = 1:numel(xx)
        t0 = xx(ii);
        if t0 > 0
            integrand = @(t) ...
                gauss_race_cdf(t0 - t, ...
                               mu(1,1), sigma(1,1), ...
                               mu(1,2), sigma(1,2))  ...
              .* gaussian_race_pdf(t, ...
                                   mu(2,1), sigma(2,1), ...
                                   mu(2,2), sigma(2,2));
            yy(ii) = integral(integrand, -Inf, Inf);
        end
    end
end


function rnd = gauss_relay_rnd(n, mu, sigma)
    rnd1 = gauss_race_rnd(n, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2));
    rnd2 = gauss_race_rnd(n, mu(2,1), sigma(2,1), mu(2,2), sigma(2,2));

    rnd = rnd1 + rnd2;
end

function millerYY = gauss_miller_cdf(xx, mu, sigma)
    F1 = gauss_uni_cdf(xx, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));
    F2 = gauss_uni_cdf(xx, mu(1,2), sigma(1,2), mu(2,2), sigma(2,2));

    millerYY = F1 + F2;
end

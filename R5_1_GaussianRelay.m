%% GETTING EXGAUSSIAN RSE + VIOLATIONS

clear; clc; close all

addpath(genpath(fullfile('Functions')));

%% Parameters

mu = [0.1 0.1;...
      0.4 0.4];
sigma = [0.05 0.05;...
         0.1  0.1];

n = 100000;
simN = 10000;
xMin = 0;
xMax = 1.1;
xx = linspace(xMin, xMax, n);

%% Test unisensory

yy = gauss_uni_cdf(xx, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));
rand = gauss_uni_rnd(n, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));

f1 = figure();
tiledlayout(1,3,'TileSpacing','compact')

nexttile
hold on

plot(xx,yy,'DisplayName','MyFunction','LineWidth',2)
plot(sort(rand),getCP(n),'LineStyle','--','DisplayName','NotMyFunction','LineWidth',2)
legend('Box','off')
title('Gaussian Unisensory')
ylabel('CDF')
xlabel('RT (s)')

%% Test gaussian race

yy = gauss_race_cdf(xx, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2));
rr = gauss_race_rnd(n, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2));

nexttile
hold on

plot(xx,yy,'DisplayName','CDF','LineWidth',2)
plot(sort(rr),getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('Gaussian Race')
ylabel('CDF')
xlabel('RT (s)')

%% Test gaussian relay

yy = gauss_relay_cdf(xx, mu, sigma);
rr = gauss_relay_rnd(n, mu, sigma);

nexttile
hold on

plot(xx,yy,'DisplayName','CDF','LineWidth',2)
plot(sort(rr),getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('Gaussian Relay')
ylabel('CDF')
xlabel('RT (s)')

%% Plot CDF of uni + multi conditions compared to RMI

aYY = gauss_uni_cdf(xx, mu(1,1), sigma(1,1), mu(2,1), sigma(2,1));
vYY = gauss_uni_cdf(xx, mu(1,2), sigma(1,2), mu(2,2), sigma(2,2));

relayYY = gauss_relay_cdf(xx, mu, sigma);
simXX = gauss_relay_rnd(simN, mu, sigma);
simYY = getCP(simN);

millerYY = gauss_miller_cdf(xx, mu, sigma);
griceYY = getGriceCDF(aYY, vYY);
violation = getViolationFromCDF(xx, relayYY, millerYY);
gain = getGainFromCDF(xx, relayYY, griceYY);

f2 = figure();
tiledlayout(1,3)
lw = 2;

nexttile
hold on

plot(xx, aYY, 'LineWidth',lw, 'DisplayName','A')
plot(xx, vYY, 'LineWidth',lw, 'DisplayName','V')
plot(sort(simXX), simYY, 'LineWidth',lw, 'DisplayName','simAV')
plot(xx, relayYY, 'LineWidth',lw, 'DisplayName','AV')
plot(xx, millerYY, 'LineWidth',lw, 'DisplayName','RMI')
xlim([xMin, xMax])
ylim([0 1])
legend('Box','off', 'Location','southeast')
ylabel('CDF')
xlabel('RT (s)')

%% Get RSE and Violations across different parameter values

muN = 20;
mus = linspace(0.2, 0.4, muN);
sigmas = linspace(0.02, 0.08, muN);

% Pre-allocate outputs
violations = zeros(1, muN);
gains      = zeros(1, muN);

% Waitbar
hWait = waitbar(0, 'Starting parallel loop...');
D = parallel.pool.DataQueue;
afterEach(D, @(val) localIncrementWaitbar(val, hWait, muN));

parfor i = 1:muN
    % muLocal = [mus(i)  mus(i);...
    %            mu(2,1) mu(2,2)];
    sigmaLocal = [sigmas(i), sigmas(i);...
                  sigma(2,1), sigma(2,2)];
    
    aYY      = gauss_uni_cdf(xx, mu(1,1), sigmaLocal(1,1), mu(2,1), sigmaLocal(2,1));
    vYY      = gauss_uni_cdf(xx, mu(1,2), sigmaLocal(1,2), mu(2,2), sigmaLocal(2,2));
    relayYY  = gauss_relay_cdf(xx, mu, sigmaLocal);
    millerYY = gauss_miller_cdf(xx, mu, sigmaLocal);
    griceYY  = getGriceCDF(aYY, vYY);

    violations(i) = getViolationFromCDF(xx, relayYY, millerYY);
    gains(i)      = getGainFromCDF(xx, relayYY, griceYY);

    % After each iteration, send '1' to indicate one more iteration is done
    send(D, 1);
end

disp('Parfor loop complete!');

nexttile
plot(mus, gains)
ylabel('RSE (s)')
xlabel('First Mu & Sigma')
nexttile
plot(mus, violations)
ylabel('Violation (s)')
xlabel('First Mu & Sigma')

%% Functions

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
    parfor ii=1:length(xx)
        if xx(ii)>0
            fun = @(t) gauss_race_cdf(xx(ii)-t, mu(1,1), sigma(1,1), mu(1,2), sigma(1,2)) .*...
                       gaussian_race_pdf(t, mu(2,1), sigma(2,1), mu(2,2), sigma(2,2));

            yy(ii) = integral( fun, -Inf, Inf );
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

% Upadte waitbar
function localIncrementWaitbar(increment, hWait, totalN)
% This function runs on the *client* side when we do send(D, increment).
% 
% We store the iteration count in a persistent variable so each call
% can "remember" how many iterations have already finished. We also
% need the total iteration count (totalN) to compute the fraction done.

    persistent count
    if isempty(count)
        count = 0;    % Initialize the first time
    end

    % Increment the count by however many iterations were just finished
    count = count + increment;

    % Fraction done
    fractionDone = count / totalN;
    
    % Update waitbar
    waitbar(fractionDone, hWait, ...
            sprintf('Processing %d of %d...', count, totalN));
    drawnow;  % Force immediate graphical update

    % If we are at or beyond total iterations, close the waitbar
    if count >= totalN
        close(hWait);
        count = 0;  % reset so next run doesn't pick up old value
    end
end

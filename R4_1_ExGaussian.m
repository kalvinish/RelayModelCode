%% GETTING EXGAUSSIAN RSE + VIOLATIONS

clear; clc; close all

addpath(genpath(fullfile('Functions')));

%% Parameters

mu = [0.4 0.4];
sigma = [0.1 0.1];
tau = [0.1 0.1];

n = 100000;
simN = 10000;
xMin = 0;
xMax = 1.1;
xx = linspace(xMin, xMax, n);

%% Test unisensory

yy = exgauss_cdf(xx, mu(1), sigma(1), tau(1));
rand = exgauss_rnd(mu(1), sigma(1), tau(1), n, 1);

f1 = figure();
tiledlayout(2,2,'TileSpacing','compact')

nexttile
hold on

plot(xx,yy,'DisplayName','MyFunction','LineWidth',2)
plot(sort(rand),getCP(n),'LineStyle','--','DisplayName','NotMyFunction','LineWidth',2)
legend('Box','off')
title('ExGaussian Unisensory')
ylabel('CDF')
xlabel('RT (s)')

%% Test exponential race

yy = exponential_race_cdf(xx, tau);
rr = exponential_race_rnd(tau, n);

nexttile
hold on

plot(xx,yy,'DisplayName','CDF','LineWidth',2)
plot(sort(rr),getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('Exponential Race')
ylabel('CDF')
xlabel('RT (s)')

%% Test gaussian race

yy = gaussian_race_pdf(xx, mu, sigma);
rr = gaussian_race_rnd(mu, sigma, n);

nexttile
hold on

plot(xx,yy,'DisplayName','PDF','LineWidth',2)
histogram(rr,100,'Normalization','pdf')
legend('Box','off')
title('Gaussian Race')
ylabel('PDF')
xlabel('RT (s)')

%% Test exgaussian relay

yy = exgauss_relay_cdf(xx, mu, sigma, tau);
rr = exgauss_relay_rnd(mu, sigma, tau, n, 1);

nexttile
hold on

plot(xx,yy,'DisplayName','CDF','LineWidth',2)
plot(rr,getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('ExGaussian Relay')
ylabel('CDF')
xlabel('RT (s)')

%% Plot CDF of uni + multi conditions compared to RMI

aYY = exgauss_cdf(xx, mu(1), sigma(1), tau(1));
vYY = exgauss_cdf(xx, mu(2), sigma(2), tau(2));

relayYY = exgauss_relay_cdf(xx, mu, sigma, tau);
simXX = exgauss_relay_rnd(mu, sigma, tau, simN, 1);
simYY = getCP(simN);

millerYY = exgauss_miller_cdf(xx, mu, sigma, tau);
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
plot(simXX, simYY, 'LineWidth',lw, 'DisplayName','simAV')
plot(xx, relayYY, 'LineWidth',lw, 'DisplayName','AV')
plot(xx, millerYY, 'LineWidth',lw, 'DisplayName','RMI')
xlim([xMin, xMax])
ylim([0 1])
legend('Box','off', 'Location','southeast')
ylabel('CDF')
xlabel('RT (s)')

%% Get RSE and Violations across different parameter values

tauN = 50;
taus = linspace(0.01, 0.4, tauN);

% Pre-allocate outputs
violations = zeros(1, tauN);
gains      = zeros(1, tauN);

% Waitbar
hWait = waitbar(0, 'Starting parallel loop...');
D = parallel.pool.DataQueue;
afterEach(D, @(val) localIncrementWaitbar(val, hWait, tauN));

parfor i = 1:tauN
    tauLocal = repelem(taus(i), 2);
    
    aYY      = exgauss_cdf(xx, mu(1), sigma(1), tauLocal(1));
    vYY      = exgauss_cdf(xx, mu(2), sigma(2), tauLocal(2));
    relayYY  = exgauss_relay_cdf(xx, mu, sigma, tauLocal);
    millerYY = exgauss_miller_cdf(xx, mu, sigma, tauLocal);
    griceYY  = getGriceCDF(aYY, vYY);

    violations(i) = getViolationFromCDF(xx, relayYY, millerYY);
    gains(i)      = getGainFromCDF(xx, relayYY, griceYY);

    % After each iteration, send '1' to indicate one more iteration is done
    send(D, 1);
end

disp('Parfor loop complete!');

nexttile
plot(taus, gains)
ylabel('RSE (s)')
xlabel('Tau')
nexttile
plot(taus, violations)
ylabel('Violation (s)')
xlabel('Tau')

%% Functions

% Race between two exponential distributions (CDF)
function yy = exponential_race_cdf(xx, tau)
    F1 = expcdf(xx, tau(1));
    F2 = expcdf(xx, tau(2));

    yy = F1 + F2 - (F1 .* F2);
end

% Race between two exponential distributions (RND)
function rr = exponential_race_rnd(tau, n)
    rr1 = exprnd(tau(1), n, 1);
    rr2 = exprnd(tau(2), n, 1);

    rr = min([rr1 rr2], [], 2);
end

% Race between two Gaussian distributions (PDF)
function yy = gaussian_race_pdf(xx, mu, sigma)
    F1 = normcdf(xx, mu(1), sigma(1));
    F2 = normcdf(xx, mu(2), sigma(2));
    f1 = normpdf(xx, mu(1), sigma(1));
    f2 = normpdf(xx, mu(2), sigma(2));
    
    yy = f1 .* (1 - F2) + f2 .* (1 - F1);
    
    yy(isnan(yy)) = 0;
end

% Race between two Gaussian distributions (RND)
function rr = gaussian_race_rnd(mu, sigma, n)
    rr1 = normrnd(mu(1), sigma(1), n, 1);
    rr2 = normrnd(mu(2), sigma(2), n, 1);

    rr = min([rr1 rr2], [], 2);
end

% Relay between exgaussian distributions (RND)
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



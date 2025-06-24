%% GETTING EXWALD RSE + VIOLATIONS

clear; clc; close all

addpath(genpath(fullfile('Functions')));

%% Parameters

mu = [0.4 0.4];
lambda = [4 4];
tau = [0.1 0.1];

n = 100000;
simN = 10000;
xMin = 0;
xMax = 1.5;
xx = linspace(xMin, xMax, n);

%% Test unisensory

tic
yy = exwald_cdf(xx, mu(1), lambda(1), tau(1));
toc

rnd = exwald_rnd(n, mu(1), lambda(1), tau(1));

f1 = figure();
tiledlayout(2,2,'TileSpacing','compact')

nexttile
hold on

plot(xx,yy,'DisplayName','MyFunction','LineWidth',2)
plot(sort(rnd),getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('ExWald Unisensory')
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

%% Test wald race

yy = wald_race_pdf(xx, mu, lambda);
rr = wald_race_rnd(mu, lambda, n);

nexttile
hold on

plot(xx,yy,'DisplayName','PDF','LineWidth',2)
histogram(rr,100,'Normalization','pdf','DisplayName','RND')
legend('Box','off')
title('Wald Race')
ylabel('PDF')
xlabel('RT (s)')

%% Test exwald relay

yy = exwald_relay_cdf(xx, mu, lambda, tau);
rr = exwald_relay_rnd(mu, lambda, tau, n, 1);

nexttile
hold on

plot(xx,yy,'DisplayName','CDF','LineWidth',2)
plot(rr,getCP(n),'LineStyle','--','DisplayName','RND','LineWidth',2)
legend('Box','off')
title('ExWald Relay')
ylabel('CDF')
xlabel('RT (s)')

%% Plot CDF of uni + multi conditions compared to RMI

aYY = exwald_cdf(xx, mu(1), lambda(1), tau(1));
vYY = exwald_cdf(xx, mu(2), lambda(2), tau(2));

relayYY = exwald_relay_cdf(xx, mu, lambda, tau);
simXX = exwald_relay_rnd(mu, lambda, tau, simN, 1);
simYY = getCP(simN);

millerYY = exwald_miller_cdf(xx, mu, lambda, tau);
griceYY = getGriceCDF(aYY, vYY);
violation = getViolationFromCDF(xx, relayYY', millerYY);
gain = getGainFromCDF(xx, relayYY', griceYY);

f2 = figure();
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

%% SAVE FIGURE
exportgraphics(gcf, fullfile(pwd, 'Figures', 'ExWaldExample.pdf'), 'ContentType','vector');

%% Functions

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


function z = exwald_rnd(n, mu, lambda, t)
%EXWALD_RND  Simulate random draws from the ex-Wald distribution
%
%   Z = EXWALD_RND(N, MU, LAMBDA, T) returns an N-by-1 vector Z of
%   independent draws from the ex-Wald distribution whose Wald (IG)
%   component has mean MU and shape LAMBDA, and whose exponential
%   component has mean T.  It simulates Z = W + X by
%      W = random('InverseGaussian',MU,LAMBDA)
%      X = exprnd(T)
%
%   Example:
%     % 10,000 draws with mu=2, lambda=3, t=0.5
%     z = exwald_rnd(10000, 2, 3, 0.5);
%     histogram(z,100);  % empirical distribution

    % Input validation
    if ~isscalar(n) || n<=0 || floor(n)~=n
        error('N must be a positive integer scalar.');
    end

    % Preallocate
    z = zeros(n,1);

    % Create an Inverse Gaussian distribution object
    pd = makedist('InverseGaussian','mu',mu,'lambda',lambda);

    % Draw
    W = random(pd, n, 1);    % Wald / IG variates
    X = exprnd(t, n, 1);     % Exponential variates with mean t

    % Sum to get ex-Wald
    z = W + X;
end

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

% Race between two Gaussian distributions (RND)
function rr = wald_race_rnd(mu, lambda, n)
    rr1 = random('inversegaussian', mu(1), lambda(1), [n,1]);
    rr2 = random('inversegaussian', mu(2), lambda(2), [n,1]);

    rr = min([rr1 rr2], [], 2);
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

% Relay between exgaussian distributions (RND)
function rr = exwald_relay_rnd(mu, lambda, tau, n, sort_flg)
    rr1 = exponential_race_rnd(tau, n);
    rr2 = wald_race_rnd(mu, lambda, n);

    if ~exist('sort_flg','var')
        rr = rr1 + rr2;
    elseif sort_flg == 1
        rr = sort(rr1 + rr2);
    else
        rr = rr1 + rr2;
    end
end

% Race Model Inequality for ExGaussian
function yy = exwald_miller_cdf(xx, mu, lambda, tau)
    F1 = exwald_cdf(xx, mu(1), lambda(1), tau(1));
    F2 = exwald_cdf(xx, mu(2), lambda(2), tau(2));
    
    yy = F1 + F2;
end


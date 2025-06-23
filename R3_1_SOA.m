% SIMULATE SOA DATA USING EMPIRICAL DATA FROM MILLER 1986 (GET WEIGHT USING RMSE)

% Clear workspace, command window, and close all figures
clear
clc
close all

addpath(fullfile(cd, "Functions"))

%% PARAMETERS FOR SIMULATION

% Participant BD data
% Fit Inverse Gaussian distribution parameters from descriptive statistics
[muA, lambdaA, res1] = fitIGfromDesc(0.231, 0.0028, 0.219, 400); % Auditory
[muV, lambdaV, res2] = fitIGfromDesc(0.348, 0.0046, 0.326, 400); % Visual

% Uncomment below for participant KY data
% [muA, lambdaA, res1] = fitIGfromDesc(0.211, 0.003, 0.193, 400);
% [muV, lambdaV, res2] = fitIGfromDesc(0.282, 0.0031, 0.266, 400);

%% SIMULATE SOA EXPERIMENT

% Number of trials
n = 400;

% SOA lags from Miller's experiment (in seconds)
millerLags = [-0.167 -0.133 -0.100 -0.067 -0.033 0 0.033 0.067 0.100 0.133 0.167];

% Simulate lags for the experiment (using 100 points between -0.167 and 0.167)
lags = linspace(-0.167, 0.167, 100);
lagsN = length(lags);
lagLabels = lags * 1000; % Convert to milliseconds

% Empirical data from Miller (1986) for participant BD
millerMeans = [234 230 227 228 221 217 238 263 277 298 316]./1000;
millerMedians = [223 222 221 223 219 215 233 259 274 296 316]./1000;
millerMeanSE = [2.9 2.0 2.0 1.6 1.4 1.4 1.4 1.3 1.5 1.6 1.7];
millerRSE = [-3 1 4 3 10 14 26 35 54 50 32]./1000;

% Uncomment below for participant KY data
% millerMeans = [216 217 214 218 215 208 237 249 256 273 278];
% millerMedians = [196 194 192 192 198 190 216 235 244 264 268];
% millerMeanSE = [3.7 3.8 3.9 3.8 3.3 3.2 3.1 2.9 2.3 2.7 3.1];
% millerRSE = [-5 -6 -3 -7 -4 3 7 29 26 9 4];

%% FIND OPTIMAL WEIGHTS USING RMSE

% Set the sample size for optimization
optim_n = 100000;

parfor i = 1:10
    % Optimize the weight parameter w in [0, 0.5] to minimize RMSE
    [w_opt_raw(i), fval(i)] = fminbnd(@(w) computeRMSE(w, optim_n, millerLags, muA, muV, lambdaA, lambdaV, millerRSE), 0, 0.5);
    disp(i)
end

% [w_opt_raw(i), fval] = fminbnd(@(w) computeRMSE(w, optim_n, muA, muV, lambdaA, lambdaV, millerRSE(6)), 0, 0.5);


w_opt = mean(w_opt_raw);

disp(w_opt)
disp(std(w_opt_raw))

% Use the optimized weight for both w1_1 and w2_1 (equal values)
w1_1 = w_opt; w1_2 = 1 - w1_1;
w2_1 = w_opt; w2_2 = 1 - w2_1;

%% RUN SIMULATION WITH OPTIMIZED WEIGHTS

% Number of repetitions for the simulation
repN = 10000;

means = nan(lagsN, repN, 3);
SEs = nan(lagsN, repN);
medians = nan(lagsN, repN, 3);
RSEs = nan(lagsN, repN);
minRT = nan(lagsN, repN);
meanRace = nan(lagsN, repN);
rseRace = nan(lagsN, repN);

% Preallocate rtA, rtV, rtAV, rtRace
rtA = nan(n, lagsN);
rtV = nan(n, lagsN);
rtAV = nan(n, lagsN);
rtRace = nan(n, lagsN);

% Initialize the waitbar to monitor progress
hWait = waitbar(0, 'Starting...');

for rep = 1:repN
    parfor i = 1:lagsN
        lag = lags(i);

        % Determine lags for auditory and visual stimuli
        if lag < 0
            lagA = 0;
            lagV = abs(lag);
        else
            lagA = abs(lag);
            lagV = 0;
        end

        % Generate reaction times for auditory and visual stimuli with lags
        rtA(:, i) = uniRND(n, muA, lambdaA, 1) + lagA;
        rtV(:, i) = uniRND(n, muV, lambdaV, 1) + lagV;

        % Generate reaction times for audiovisual stimuli using relay model
        rtAV(:, i) = lagRND(n, muA, muV, lambdaA, lambdaV, w1_1, w1_2, w2_1, w2_2, lag, 1);

        % Generate reaction times for audiovisual stimuli using race model
        rtRace(:, i) = lagRND(n, muA, muV, lambdaA, lambdaV, 0, 1, 0, 1, lag, 1);

        % Compute minimum mean RT between auditory and visual
        minRT(i, rep) = min([mean(rtA(:, i)), mean(rtV(:, i))]);

        % Compute RSE (Redundancy Gain)
        RSEs(i, rep) = (minRT(i, rep) - mean(rtAV(:, i))) * 1000;

        % Store statistics
        % means(i, rep) = mean(rtAV(:, i)) * 1000;
        SEs(i, rep) = standard_error(rtAV(:, i)) * 1000;

        % Assign medians in a single statement
        medians(i, rep, :) = [median(rtA(:, i)), median(rtV(:, i)), median(rtAV(:, i))] * 1000;
        means(i, rep, :) = [mean(rtA(:, i)), mean(rtV(:, i)), mean(rtAV(:, i))] * 1000;

        meanRace(i, rep) = mean(rtRace(:, i)) * 1000;
        rseRace(i, rep) = (minRT(i, rep) - mean(rtRace(:, i))) * 1000;
    end

    % Update the waitbar after each iteration
    waitbar(rep / repN, hWait, sprintf('Processing %d of %d...', rep, repN));
end

% Close the waitbar after the loop is finished
close(hWait);

%% CLEAN SUMMARY DATA

% Compute 95% confidence intervals
RSEs95 = computeQuantileCI(RSEs, 2);
rseRace95 = computeQuantileCI(rseRace, 2);

medians95 = nan(2, lagsN, 3);
means95 = nan(2, lagsN, 3);

for i = 1:3
    medians95(:,:,i) = computeQuantileCI(medians(:,:,i), 2);
    means95(:, :, i) = computeQuantileCI(means(:,:,i), 2);
end

% Transpose for plotting
RSEs95 = RSEs95';
rseRace95 = rseRace95';

% Compute overall means and standard errors
means = mean(means, 2);
SEs = mean(SEs, 2);
medians = mean(medians, 2);
RSEs = mean(RSEs, 2);
meanRace = mean(meanRace, 2);
rseRace = mean(rseRace, 2);

%% Mean Response Times

f1 = figure;
t = tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'none', 'Units', 'centimeters', 'OuterPosition', [0 0 14 5.5]);

aud_col = '#009FE3';
vis_col = '#3AAA35';
relay_col = '#706F6F';
my_linewidth = 1.5;
my_fontsize = 9;
my_markersize = 15;

ax1 = nexttile(t);
% Create shaded area for simulated median 95% CI
x2 = [lagLabels, fliplr(lagLabels)];
inBetween = [means95(1,:,3), fliplr(means95(2,:,3))];
fill(x2, inBetween, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% CI');
hold on;

lagsA = zeros(length(lags), 1);
lagsA(lags > 0) = lags(lags > 0);
aData = (muA + lagsA) * 1000;

lagsV = zeros(length(lags), 1);
lagsV(lags < 0) = abs(lags(lags < 0));
vData = (muV + lagsV) * 1000;

plot(lagLabels, aData, 'Color', aud_col, 'LineStyle', '--', "LineWidth", my_linewidth, 'DisplayName', 'A')
plot(lagLabels, vData, 'Color', vis_col, 'LineStyle', '--', "LineWidth", my_linewidth, 'DisplayName', 'V')
plot(lagLabels, means(:,:,3), "-k", "LineWidth", my_linewidth, 'DisplayName', 'AV')

plot(lagLabels, meanRace, "--r", "LineWidth", my_linewidth, 'DisplayName', 'Race')

scatter(millerLags*1000, millerMeans*1000, my_markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', my_linewidth, 'DisplayName', 'Miller (1986)')

ylim([200 400])

xlabel('SOA (ms)');
ylabel('Response Time (ms)');
legend('Location', 'northwest', 'Box', 'Off');
set(gca, 'TickDir', 'out');
box off;
ax = gca;
ax.TickDir = "Out";
ax.YColor = "k";
ax.XColor = "k";
ax.LineWidth = my_linewidth;
ax.FontSize = my_fontsize;
yticks(linspace(200, 400, 5))
xticks(linspace(-200, 200, 5))
hold off;

%% PLOT SIMULATED RSE WITH EMPIRICAL DATA

ax2 = nexttile(t);

% Create shaded area for simulated RSEs 95% CI
x2 = [lagLabels, fliplr(lagLabels)];
inBetween = [RSEs95(:,1)', fliplr(RSEs95(:,2)')];
fill(x2, inBetween, '', 'FaceColor', relay_col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;

% Create shaded area for simulated RSEs 95% CI
x2 = [lagLabels, fliplr(lagLabels)];
inBetween = [rseRace95(:,1)', fliplr(rseRace95(:,2)')];
fill(x2, inBetween, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot simulated relay model RSE
plot(lagLabels, RSEs, 'Color', relay_col, 'LineWidth', my_linewidth);

% Plot Miller's empirical RSE
scatter(millerLags*1000, millerRSE.*1000, my_markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', my_linewidth, 'DisplayName', 'Miller (1986)')

% Plot simulated race model RSE
plot(lagLabels, rseRace, 'Color', 'r', 'LineStyle', '-', 'LineWidth', my_linewidth);

xlabel('SOA (ms)');
ylabel('RSE (ms)');
% title(['w_{opt} = ' num2str(w1_1, '%.4f')]);
set(gca, 'TickDir', 'out');
yline(0, '--k', 'LineWidth', 2);
box off;
ax = gca;
ax.TickDir = "Out";
ax.YColor = "k";
ax.XColor = "k";
ax.LineWidth = my_linewidth;
ax.FontSize = my_fontsize;
hold off;

%% SAVE PLOT

exportgraphics( f1, fullfile( cd, "Figures", "figure4_fit2zero.pdf" ), "ContentType", "vector" )

%% FUNCTIONS

% Function to compute standard error
function se = standard_error(x, dim)
    % Calculates the standard error along the specified dimension
    if nargin < 2
        dim = find(size(x) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end
    std_x = std(x, 0, dim);
    n = size(x, dim);
    se = std_x / sqrt(n);
end

% Function to compute RMSE between simulated data and empirical data
function rmse = computeRMSE(w, n, lags, muA, muV, lambdaA, lambdaV, millerRSE)
    lagsN = length(lags);

    % Set weights based on input w
    if length(w) > 1
        w1_1 = w(1); w1_2 = 1 - w1_1;
        w2_1 = w(2); w2_2 = 1 - w2_1;
    else
        w1_1 = w; w1_2 = 1 - w1_1;
        w2_1 = w; w2_2 = 1 - w2_1;
    end

    for i = 1:lagsN
        lag = lags(i);

        % Determine lags for auditory and visual stimuli
        if lag < 0
            lagA = 0;
            lagV = abs(lag);
        else
            lagA = abs(lag);
            lagV = 0;
        end

        % Generate reaction times for auditory and visual stimuli with lags
        rtA = uniRND(n, muA, lambdaA, 1) + lagA;
        rtV = uniRND(n, muV, lambdaV, 1) + lagV;

        % Generate reaction times for audiovisual stimuli using relay model
        rtAV = lagRND(n, muA, muV, lambdaA, lambdaV, w1_1, w1_2, w2_1, w2_2, lag, 1);

        % Compute minimum mean RT between auditory and visual
        minRT = min([mean(rtA), mean(rtV)]);

        % Generate reaction times for audiovisual stimuli using relay model
        RSEs(i) = (minRT - mean(rtAV));
    end

    % Combine differences in means and medians
    differences = (millerRSE - RSEs);

    % Calculate RMSE
    rmse = sqrt(mean(differences .^ 2));
end

% % Function to compute RMSE between simulated data and empirical data
% function rmse = computeRMSE(w, n, muA, muV, lambdaA, lambdaV, millerRSE)
% 
%     % Set weights based on input w
%     if length(w) > 1
%         w1_1 = w(1); w1_2 = 1 - w1_1;
%         w2_1 = w(2); w2_2 = 1 - w2_1;
%     else
%         w1_1 = w; w1_2 = 1 - w1_1;
%         w2_1 = w; w2_2 = 1 - w2_1;
%     end
% 
%     % Generate reaction times for auditory and visual stimuli with lags
%     rtA = uniRND(n, muA, lambdaA, 1);
%     rtV = uniRND(n, muV, lambdaV, 1);
% 
%     % Generate reaction times for audiovisual stimuli using relay model
%     % rtAV = lagRND(n, muA, muV, lambdaA, lambdaV, w1_1, w1_2, w2_1, w2_2, 0, 1);
%     rtAV = multiRND(n, muA, muV, lambdaA, lambdaV, w1_1, w1_2, w2_1, w2_2);
% 
%     % Compute minimum mean RT between auditory and visual
%     minRT = min([mean(rtA), mean(rtV)]);
% 
%     % Generate reaction times for audiovisual stimuli using relay model
%     RSE = (minRT - mean(rtAV)) * 1000;
% 
%     % Combine differences in means and medians
%     differences = (millerRSE - RSE);
% 
%     % Calculate RMSE
%     rmse = sqrt(mean(differences .^ 2));
% end

% Function to compute 95% confidence intervals using quantiles
function ci = computeQuantileCI(data, dim)
    % Computes the 2.5% and 97.5% quantiles along the specified dimension
    %
    % Inputs:
    %   data - Numeric matrix
    %   dim  - Dimension along which to compute quantiles (1 for columns, 2 for rows)
    %
    % Output:
    %   ci - Matrix containing the lower and upper quantiles along the specified dimension

    % Validate inputs
    if nargin < 2
        error('You must specify both the data matrix and the dimension (1 for columns, 2 for rows).');
    end

    % Check if data is numeric
    if ~isnumeric(data)
        error('Input data must be a numeric matrix.');
    end

    % Check if dim is valid
    if dim ~= 1 && dim ~= 2
        error('Dimension dim must be either 1 (columns) or 2 (rows).');
    end

    % Define the quantiles to compute
    quantiles = [0.025, 0.975];

    % Compute the quantiles along the specified dimension
    ci_lower = quantile(data, quantiles(1), dim);
    ci_upper = quantile(data, quantiles(2), dim);

    % Concatenate the lower and upper quantiles
    if dim == 1
        ci = [ci_lower, ci_upper]';
    else
        ci = [ci_lower'; ci_upper'];
    end
end


%% Figure 2: Increasing Number of Stages

clear
clc
close all

addpath(genpath(fullfile('Functions')));

%% Set Parameters

load(fullfile(cd, "FittedParameters", "Miller82_Parameters.mat"))

A = params{1,1};
V = params{1,2};
W = params{1,3};

aMU = A.mu;
vMU = V.mu;
aLAMBDA = A.lambda;
vLAMBDA = V.lambda;

xMin = 0;
xMax = 700;
xN = 1000;

xx = linspace(xMin, xMax, xN);

raceMax = 10;
raceN = raceMax;

%% LOAD DATA TO WORKSPACE

loadPath = fullfile(cd, "EmpiricalData", "Miller82");

data = nan(10, 3);
params = nan(3,2);

for i = 1:3
    temp_data = readmatrix(fullfile(loadPath, [int2str(i), '.csv']));
    data(:,i) = temp_data(:,1);
end

%% Get CDFs

multi = nan(xN, raceN);

races = linspace(1,raceMax,raceN);

parfor i = 1:raceN
    race = races(i);
    multi(:,i) = multipleRacesCDF(xx, race, aMU, vMU, aLAMBDA, vLAMBDA);
    disp(i)
end

%% Plot CDFs

f1 = figure;
t = tiledlayout(1, 4, 'TileSpacing', 'none', 'Padding', 'none', 'Units', 'centimeters', 'OuterPosition', [0 0 16 5.5]);
ax1 = nexttile(t, [1, 2]);

my_linewidth = 1.5;
my_fontsize = 9;

hold on;

my_linewidth = 1.5;
my_fontsize = 9;

aud_col = '#009FE3';
vis_col = '#3AAA35';

uniA = uniCDF(xx, aMU, aLAMBDA);
uniV = uniCDF(xx, vMU, vLAMBDA);

plot(xx, uniA, 'Color', aud_col, 'LineWidth', my_linewidth)
plot(xx, uniV, 'Color', vis_col, 'LineWidth', my_linewidth)

races = 1:raceN;

% Define a colormap using cbrewer2
colors = cbrewer2('Reds', length(races));

% Normalize weights to be within the range of the colormap
norm_weights = (races - min(races)) / (max(races) - min(races));
color_indices = round(norm_weights * (size(colors, 1) - 1)) + 1;

colorMapLength = raceN;
red = [0.8 0.8 0.8];
pink = [0 0 0];
colors_p = [linspace(red(1),pink(1),colorMapLength)', linspace(red(2),pink(2),colorMapLength)', linspace(red(3),pink(3),colorMapLength)'];

% Plot each CDF with the corresponding color
for i = 1:raceN
    plot(xx, multi(:,i), 'Color', colors_p(color_indices(i), :), 'LineWidth', my_linewidth);
end

n = 100000;
miller = getMiller([uniRND(n, aMU, aLAMBDA) uniRND(n, vMU, vLAMBDA)]);
raab = getRaab([uniRND(n, aMU, aLAMBDA) uniRND(n, vMU, vLAMBDA)]);
plot(miller, getCP(n), 'r', 'LineWidth', 2, 'LineStyle', '-')
plot(raab, getCP(n), 'r', 'LineWidth', 2, 'LineStyle', '--')

dataCDFvals = linspace(0.05, 0.95, 10);

scatter(data(:,1), dataCDFvals, 15, 'MarkerEdgeColor', aud_col, 'MarkerFaceColor', 'w', 'LineWidth', 1.5)
scatter(data(:,2), dataCDFvals, 15, 'MarkerEdgeColor', vis_col, 'MarkerFaceColor', 'w', 'LineWidth', 1.5)
scatter(data(:,3), dataCDFvals, 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5)

% Create a custom colorbar
colormap(colors_p); % Use the defined colors
c = colorbar('eastoutside', 'Visible', 'on');
clim([min(races) max(races)]);
ylabel(c, "Race Number")

% Add labels and title
xlabel('Response Time (ms)');
ylabel('Cumulative Probability');
hold off;
set(gca,'TickDir','out');
set(gca, 'XColor', 'k');
set(gca, 'YColor', 'k');
xlim([100 xMax])
ylim([0 1])
set(gca, 'linewidth', my_linewidth)
set(gca, 'FontSize', my_fontsize)
yticks(linspace(0, 1, 3))
xticks(linspace(100, 700, 3))

% Get percent increase in RSE with increasing race number

millerCDFvals = millerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)';
raabCDFvals = raabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)';
griceCDF = getGriceCDF(uniA, uniV)';

for i = 1:raceN
    rse(i) = getGainFromCDF(xx, multi(:,i), griceCDF);
    violation(i) = getViolationFromCDF(xx, multi(:,i), millerCDFvals);
end

diffRSE = diff(rse);

pct_increase = percentageIncrease(rse)

% PLOTS

% Create the plot
ax2 = nexttile(t);
hold on
plot(1, rse(1), 'LineWidth', my_linewidth, 'Marker', 'o', 'MarkerFaceColor', 'w', 'Color', 'r', 'MarkerSize', 4)
plot(2:raceN, rse(2:end), 'LineWidth', my_linewidth, 'Marker', 'o', 'MarkerFaceColor', 'w', 'Color', 'k', 'MarkerSize', 4)
ylabel('RSE (ms)')
xlabel('Races (#)')
set(gca, 'TickDir', 'out');
set(gca, 'linewidth', my_linewidth)
set(gca, 'FontSize', my_fontsize)
set(gca, 'XColor', 'k');
set(gca, 'YColor', 'k');
box off
xlim([0 raceMax])
xticks(linspace(0, raceMax, raceN/2))
ylim([0 200])
yline(rse(1), 'LineStyle', '--', 'Color', 'r', 'LineWidth', my_linewidth)

ax3 = nexttile(t);
hold on
plot(1, violation(1), 'LineWidth', my_linewidth, 'Marker', 'o', 'MarkerFaceColor', 'w', 'Color', 'r', 'MarkerSize', 4)
plot(2:raceN, violation(2:end), 'LineWidth', my_linewidth, 'Marker', 'o', 'MarkerFaceColor', 'w', 'Color', 'k', 'MarkerSize', 4)
ylabel("Violation (ms)")
xlabel("Races (#)")
set( gca, 'TickDir', 'out' );
set( gca, 'linewidth', my_linewidth)
set( gca, "FontSize", my_fontsize)
set(gca, 'XColor', 'k');
set(gca, 'YColor', 'k');
ylim([0 100])
xlim([0 raceMax])
xticks(linspace(0, raceMax, 5))
yline(violation(1), 'LineStyle', '--', 'Color', 'r', 'LineWidth', my_linewidth)
box off

%% Save Figure
exportgraphics(f1, fullfile(cd, "Figures", "figure3.pdf"), "ContentType", "vector")

%% FUNCTIONS

function pct_increase = percentageIncrease(vec)
    % Validate that the input vector has more than one element
    if length(vec) < 2
        error('Input vector must have at least two elements.');
    end
    
    % Calculate the percentage increase
    pct_increase = ((vec(2:end) - vec(1:end-1)) ./ vec(1:end-1)) * 100;
end

%% Figure 2: Increasing Weights and the Effect on Violations

clear
clc
close all

addpath(fullfile(cd, "Functions"))

load(fullfile(cd, "FittedParameters", "Miller82_Parameters.mat"))

%% Set Parameters

A = params{1,1};
V = params{1,2};
W = params{1,3};

aMU = A.mu;
vMU = V.mu;
aLAMBDA = A.lambda;
vLAMBDA = V.lambda;

optimW = W.optimal_w;

xN = 100;

xx = linspace(0.1, 0.8, xN);

N = 100;

weights = linspace(0, 0.5, N);

colorMapLength = N;
grey = [228 210 231] / 255;
black = [137 41 133] / 255;
colors_p = [linspace(grey(1),black(1),colorMapLength)', linspace(grey(2),black(2),colorMapLength)', linspace(grey(3),black(3),colorMapLength)'];

%% LOAD DATA TO WORKSPACE

loadPath = fullfile(cd, "EmpiricalData", "Miller82");

data = nan(10, 3);

for i = 1:3
    loadedData = readmatrix(fullfile(loadPath, [int2str(i), '.csv']));
    data(:,i) = loadedData(:,1);
end

dataCDFvals = linspace(0.05, 0.95, 10);

%% Plot CDF

yy = nan(xN, N);

parfor i = 1:N
    yy(:,i) = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, weights(i), 1-weights(i),...
                                                  weights(i), 1-weights(i));
end

%% Create a figure for the CDF plot

% Create your figure
f1 = figure;
t = tiledlayout(1, 4, 'TileSpacing', 'none', 'Padding', 'none', 'Units', 'centimeters', 'OuterPosition', [0 0 16 5.5]);
ax1 = nexttile(t, [1, 2]);
my_linewidth = 1.5;
my_fontsize = 9;

aud_col = '#009FE3';
vis_col = '#3AAA35';

axesHandle = ax1;
hold on;

uniA = uniCDF(xx, aMU, aLAMBDA);
uniV = uniCDF(xx, vMU, vLAMBDA);

plot(ax1, xx*1000, uniA, 'Color', aud_col, 'LineWidth', my_linewidth);
plot(ax1, xx*1000, uniV, 'Color', vis_col, 'LineWidth', my_linewidth);

% Define a colormap using cbrewer2
colors = cbrewer2('Reds', length(weights));

% Normalize weights to be within the range of the colormap
norm_weights = (weights - min(weights)) / (max(weights) - min(weights));
color_indices = round(norm_weights * (size(colors_p, 1) - 1)) + 1;

% Plot each CDF with the corresponding color
for i = 1:N
    plot(ax1, xx*1000, yy(:,i), 'Color', colors_p(color_indices(i), :), 'LineWidth', my_linewidth);
end

miller_xx = linspace(0.1, 0.8, 500);
miller = millerCDF(miller_xx, aMU, vMU, aLAMBDA, vLAMBDA);
plot(ax1, miller_xx*1000, miller, 'r', 'LineWidth', my_linewidth, 'LineStyle', '-');

raab = raabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
plot(ax1, xx*1000, raab, 'r', 'LineWidth', my_linewidth, 'LineStyle', '--');

fittedMultiCDF = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, optimW, 1-optimW, optimW, 1-optimW);
plot(ax1, xx*1000, fittedMultiCDF, 'k', 'LineWidth', my_linewidth);

scatter(ax1, data(:,1), dataCDFvals, 15, 'MarkerEdgeColor', aud_col, 'MarkerFaceColor', 'w', 'LineWidth', 1.5)
scatter(ax1, data(:,2), dataCDFvals, 15, 'MarkerEdgeColor', vis_col, 'MarkerFaceColor', 'w', 'LineWidth', 1.5)
scatter(ax1, data(:,3), dataCDFvals, 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5)

% Create a custom colorbar
colormap(colors_p); % Use the defined colors
c = colorbar('eastoutside', 'Visible','off');
clim([min(weights) max(weights)]);
% ylabel(c, "1st Race Proportion")

% Add labels and title
xlabel('Response Time (ms)');
ylabel('Cumulative Probability');
set(gca,'TickDir','out');
xlim([100 700])

% Set axis ticks and other properties
yticks(linspace(0, 1, 3))
xticks(linspace(100, 700, 3))

% ax = gca;
ax1.TickDir = "Out";
ax1.YColor = "k";
ax1.XColor = "k";
ax1.LineWidth = my_linewidth;
ax1.FontSize = my_fontsize;

%% Get RSE Plot

if 1
    parfor i = 1:N
        griceCDF = getGriceCDF(uniCDF(xx, aMU, aLAMBDA), uniCDF(xx, vMU, vLAMBDA));
        raabCDF = myRaabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
        relayCDF = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, weights(i), 1-weights(i),...
                                        weights(i), 1-weights(i));
        millercdf = millerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);

        gains(i) = getGainFromCDF(xx, relayCDF, griceCDF);
        raabGain(i) = getGainFromCDF(xx, raabCDF, griceCDF);
        violationsTrue(i) = getViolationFromCDF(xx, relayCDF, millercdf);
    end
    
    % nexttile
    % subplot(1,2,1)
    ax2 = nexttile(t);
    hold on
    % fill([weights fliplr(weights)], [conf_intervals(:,1)'*1000 fliplr(conf_intervals(:,2)'*1000)], [0.9 0.9 0.9], 'EdgeColor', 'none');
    plot(weights, gains*1000, 'k', 'LineWidth', my_linewidth)
    ylabel("RSE (ms)")
    xlabel("RT Share (%)")
    set(gca,'TickDir','out');
    box off
    xlim([0 0.5])
    ylim([0 80])
    yline(raabGain(1)*1000, '--r', 'LineWidth', my_linewidth)

    % Gains increase:
    disp(gains(end) / gains(1))
    
    % Set axis ticks and other properties
    yticks(linspace(0, 80, 5))
    xticks(linspace(0, 0.50, 6))
    xticklabels(linspace(0, 50, 6))
    
    ax2.TickDir = "Out";
    ax2.YColor = "k";
    ax2.XColor = "k";
    ax2.LineWidth = my_linewidth;
    ax2.FontSize = my_fontsize;
        
    % nexttile
    % subplot(1,2,2)
    ax3 = nexttile(t);
    hold on
    
    % Plot the main curve
    % fill([weights fliplr(weights)], [conf_intervals(:,1)'*1000 fliplr(conf_intervals(:,2)'*1000)], [0.9 0.9 0.9], 'EdgeColor', 'none');
    plot(weights, violationsTrue*1000, 'k', 'LineWidth', my_linewidth)
    ylabel("Violation (ms)")
    xlabel("RT Share (%)")
    set(gca,'TickDir','out');
    ylim([0 10])
    xlim([0 0.5])
    box off

    % Violations increase:
    disp(violationsTrue(end) / violationsTrue(2))
    
    % Set axis ticks and other properties
    yticks(linspace(0, 10, 6))
    xticks(linspace(0, 0.5, 6))
    xticklabels(linspace(0, 50, 6))
    
    yline(0, "LineStyle", "--", "Color", "r")
    
    ax3.TickDir = "Out";
    ax3.YColor = "k";
    ax3.XColor = "k";
    ax3.LineWidth = my_linewidth;
    ax3.FontSize = my_fontsize;

end

%% SAVE
exportgraphics(f1, fullfile(cd, "Figures", "figure2.pdf"), "ContentType", "vector")
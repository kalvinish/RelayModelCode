%% Function to write struct for plotting options

function plotOpts = createPlotOpts(xlim, xtickN, ylim, ytickN)

% Options for plotting
plotOpts.linewidth = 2;
plotOpts.fontsize = 14;
plotOpts.markersize = 7;
plotOpts.markerlinewidth = 2.5;
plotOpts.markerfacecol = 'w';
plotOpts.tickdir = 'out';
% Plotting colours
plotOpts.audCol = '#3AAA35';
plotOpts.visCol = '#009FE3';
plotOpts.audvisCol = '#1d1d1b';
plotOpts.modelCol = '#e52421';
plotOpts.axisCol = 'k';

plotOpts.xlim = xlim;
plotOpts.xticks = linspace(plotOpts.xlim(1),plotOpts.xlim(2),xtickN);
plotOpts.ylim = ylim;
plotOpts.yticks = linspace(plotOpts.ylim(1),plotOpts.ylim(2),ytickN);

end
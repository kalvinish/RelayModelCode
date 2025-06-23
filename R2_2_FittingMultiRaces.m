%% Fitting Miller 1982 Data with Relay Race (number of races)

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

%% LOAD DATA TO WORKSPACE

loadPath = fullfile(cd, "EmpiricalData", "Miller82");

data = nan(10, 3);

for i = 1:3
    temp_data = readmatrix(fullfile(loadPath, [int2str(i), '.csv']));
    data(:,i) = temp_data(:,1);
end

%% Get CDFs

xMin = 0;
xMax = 700;
xN = 1000;
% xx = linspace(xMin, xMax, xN);

empirical_F = linspace(0.05, 0.95, 10)';
xx = data(:,3);

raceMax = 3;
races = 1:raceMax;
raceN = length(races);

multi = nan(xN, raceN);

parfor i = 1:raceN
    race = races(i);
    [rmse(i), F(:,i)] = getRMSE(race, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA);
    disp(i)
end

%% Plot fits

f1 = figure();
hold on

scatter(data(:,3), empirical_F, 20, 'r', 'filled')

plot(xx, F(:,1))
plot(xx, F(:,2))
plot(xx, F(:,3))

full_xx = linspace(xMin, xMax, xN);
testF_1 = getMultiRelayCDF(full_xx, 1, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_1, "Color", "r", "LineStyle", ':')
testF_2 = getMultiRelayCDF(full_xx, 2, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_2, "Color", "r", "LineStyle", '--')
testF_3 = getMultiRelayCDF(full_xx, 3, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_3, "Color", "r", "LineStyle", ':')

%% Subfunction to compute the RMSE between CDFs
function [dist, predicted_F] = getRMSE(race, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA)
    % Compute predicted CDF using multiCDF
    predicted_F = getMultiRelayCDF(xx, race, aMU, vMU, aLAMBDA, vLAMBDA);

    % Compute residuals
    residuals = empirical_F - predicted_F;
    
    % Compute root mean square error (RMSE)
    dist = sqrt(mean(residuals.^2));
end

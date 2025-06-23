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
    data(:,i) = temp_data(:,1) / 1000;
end

%% Get CDFs

xMin = 0;
xMax = 0.7;
xN = 1000;
% xx = linspace(xMin, xMax, xN);

empirical_F = linspace(0.05, 0.95, 10)';
xx = data(:,3);

raceMax = 10;
races = 1:raceMax;
raceN = length(races);

multi = nan(xN, raceN);

parfor i = 1:raceN
    race = races(i);
    [rmse(i), F(:,i)] = compute_distance(race, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA);
    disp(i)
end

%% Plot fits

f1 = figure();
hold on

scatter(data(:,3), empirical_F, 20, 'r', 'filled')

plot(xx, F(:,1))
plot(xx, F(:,2))
plot(xx, F(:,3))
plot(xx, F(:,4))
plot(xx, F(:,5))

full_xx = linspace(xMin, xMax, xN);
testF_1 = multipleRacesCDF(full_xx, 1, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_1, "Color", "r", "LineStyle", '--')
testF_2 = multipleRacesCDF(full_xx, 2, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_2, "Color", "r", "LineStyle", '--')
testF_3 = multipleRacesCDF(full_xx, 3, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_3, "Color", "r", "LineStyle", '--')
testF_4 = multipleRacesCDF(full_xx, 4, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_4, "Color", "r", "LineStyle", '--')
testF_5 = multipleRacesCDF(full_xx, 5, aMU, vMU, aLAMBDA, vLAMBDA);
plot(full_xx, testF_5, "Color", "r", "LineStyle", '--')

%% Subfunction to compute the RMSE between CDFs
function [dist, predicted_F] = compute_distance(race, xx, empirical_F, aMU, vMU, aLAMBDA, vLAMBDA)
    % Compute predicted CDF using multiCDF
    predicted_F = multipleRacesCDF(xx, race, aMU, vMU, aLAMBDA, vLAMBDA);

    % Compute residuals
    residuals = empirical_F - predicted_F;
    
    % Compute root mean square error (RMSE)
    dist = sqrt(mean(residuals.^2));
end

function yy = multipleRacesCDF(xx, raceN, aMU, vMU, aLAMBDA, vLAMBDA)
    % For a single race, use raab's race.
    if raceN == 1
        yy = myRaabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
        yy = yy + randn(size(yy)) * 1e-100 + 1e-100;
        return
    end

    % Define stage weights (both stage1 and stage2)
    aW1 = 1 / raceN;
    vW1 = 1 / raceN;
    aW2 = 1 / raceN;
    vW2 = 1 / raceN;
    
    % Determine original range and extend it by a buffer.
    a_original = min(xx);
    b_original = max(xx);
    % Choose a buffer value based on the characteristic scale of the distributions.
    buffer = 5; % You may need to adjust this value.
    
    % Create an extended dense grid.
    extendedMin = max(0, a_original - buffer);
    extendedMax = b_original + buffer;
    denseX = linspace(extendedMin, extendedMax, 1000);
    
    % Compute the initial stage1 CDF on the extended dense grid.
    currentCDF_dense = stage1CDF(denseX, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1);
    
    % Loop over race stages (starting from race = 2)
    for race = 2:raceN
        newCDF_dense = zeros(size(denseX));
        numPoints = numel(denseX);
        currentStage = race;   % capture current stage
        
        % Create a waitbar if figures can be shown; otherwise, use console output.
        if usejava('desktop') && feature('ShowFigureWindows')
            hWait = waitbar(0, sprintf('Race Stage %d: 0%% complete', currentStage));
        else
            hWait = []; % No graphical display available.
            fprintf('Race Stage %d: 0%% complete\n', currentStage);
        end
        
        % Create a DataQueue for progress reporting.
        dq = parallel.pool.DataQueue();
        afterEach(dq, @(dummy) updateProgress(dummy, currentStage, numPoints, hWait));
        
        % Compute the convolution at each point in the extended grid.
        parfor ii = 1:numPoints
            if denseX(ii) > 0
                % Define an anonymous function that does the integrand:
                % It evaluates the current CDF at (denseX(ii)-t) and multiplies by the stage2PDF at t.
                fun = @(t) interp1(denseX, currentCDF_dense, denseX(ii)-t, 'linear', 0) .* ...
                           stage2PDF(t, aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2);
                % Compute the integral.
                newValue = integral(fun, 0, Inf, 'AbsTol',1e-8, 'RelTol',1e-6);
                newCDF_dense(ii) = newValue;
            else
                newCDF_dense(ii) = 0;
            end
            send(dq, 1);
        end
        
        % Close the waitbar if it was created.
        if ~isempty(hWait) && ishandle(hWait)
            close(hWait);
        end
        
        % Enforce monotonicity
        newCDF_dense = cummax(newCDF_dense);
        % Update the current CDF for the next stage.
        currentCDF_dense = newCDF_dense;
    end
    
    % Interpolate the result from the extended grid back onto the original xx values.
    yy = interp1(denseX, currentCDF_dense, xx, 'linear', 'extrap');
    yy = yy + randn(size(yy)) * 1e-100 + 1e-100;
end

function F = stage1CDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1)
    % Compute two Inverse Gaussian CDFs and combine them.
    F1 = cdf("InverseGaussian", xx, aMU * aW1, aLAMBDA * aW1^2);
    F2 = cdf("InverseGaussian", xx, vMU * vW1, vLAMBDA * vW1^2);
    F = F1 + F2 - (F1 .* F2);
    F(isnan(F)) = 0;
end

function p = stage2PDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2)
    % Compute two Inverse Gaussian PDFs and associated CDFs then combine.
    F1 = cdf("InverseGaussian", xx, aMU * aW2, aLAMBDA * aW2^2);
    F2 = cdf("InverseGaussian", xx, vMU * vW2, vLAMBDA * vW2^2);
    f1 = pdf("InverseGaussian", xx, aMU * aW2, aLAMBDA * aW2^2);
    f2 = pdf("InverseGaussian", xx, vMU * vW2, vLAMBDA * vW2^2);
    p = f1 .* (1 - F2) + f2 .* (1 - F1);
    p(isnan(p)) = 0;
end

% Local function updateProgress
function updateProgress(~, stage, numPoints, hWait)
    % Use a persistent container (Map) to hold progress for each stage.
    persistent progressCounters
    if isempty(progressCounters)
        progressCounters = containers.Map('KeyType','char','ValueType','double');
    end
    key = sprintf('%d', stage);
    if ~isKey(progressCounters, key)
        progressCounters(key) = 0;
    end
    progressCounters(key) = progressCounters(key) + 1;
    fractionDone = progressCounters(key) / numPoints;
    msg = sprintf('Race Stage %d: %.1f%% complete', stage, fractionDone * 100);
    
    if ~isempty(hWait) && ishandle(hWait)
        waitbar(fractionDone, hWait, msg);
    else
        % If no waitbar is available, output to command window.
        fprintf('%s\n', msg);
    end
    
    % Reset the counter when complete.
    if progressCounters(key) >= numPoints
        remove(progressCounters, key);
    end
end






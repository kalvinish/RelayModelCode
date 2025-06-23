function yy = getMultiRelayCDF(xx, raceN, aMU, vMU, aLAMBDA, vLAMBDA)
% getMultiRelayCDF  Compute relay-model CDF across multiple race stages
%
%   yy = getMultiRelayCDF(xx, raceN, aMU, vMU, aLAMBDA, vLAMBDA) returns the
%   cumulative distribution function (CDF) for a relay (multi-stage race)
%   model with 'raceN' sequential races. For raceN = 1, it delegates to
%   myRaabCDF (Raab's independent race model) and adds a negligible jitter.
%
%   Inputs:
%     xx       - Vector of response-time values (increasing)
%     raceN    - Number of sequential race stages (positive integer)
%     aMU      - Auditory Inverse-Gaussian mean
%     vMU      - Visual   Inverse-Gaussian mean
%     aLAMBDA  - Auditory Inverse-Gaussian lambda
%     vLAMBDA  - Visual   Inverse-Gaussian lambda
%
%   Output:
%     yy       - CDF values at xx under the relay model
%
    % Validate inputs
    validateattributes(xx,      {'numeric'}, {'vector','real','increasing'}, mfilename, 'xx', 1);
    validateattributes(raceN,   {'numeric'}, {'scalar','integer','positive'},       mfilename, 'raceN', 2);
    validateattributes(aMU,     {'numeric'}, {'scalar','real','positive'},         mfilename, 'aMU', 3);
    validateattributes(vMU,     {'numeric'}, {'scalar','real','positive'},         mfilename, 'vMU', 4);
    validateattributes(aLAMBDA, {'numeric'}, {'scalar','real','positive'},         mfilename, 'aLAMBDA', 5);
    validateattributes(vLAMBDA, {'numeric'}, {'scalar','real','positive'},         mfilename, 'vLAMBDA', 6);

    % For a single race, use Raab's independent-race CDF
    if raceN == 1
        yy = getRaabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
        yy = yy + eps * randn(size(yy));
        return;
    end

    % Stage weights: equal division across raceN stages
    w1 = 1 / raceN;
    w2 = w1;

    % Extend grid for convolution
    buffer = 200;
    xMin   = max(0, min(xx) - buffer);
    xMax   =      max(xx) + buffer;
    denseX = linspace(xMin, xMax, 10000);
    nDense = numel(denseX);

    % Initial-stage CDF
    currentCDF = stage1CDF(denseX, aMU, vMU, aLAMBDA, vLAMBDA, w1, w1);

    % Loop through additional stages
    for stage = 2:raceN
        newCDF = zeros(size(denseX));

        % Optional progress bar
        useGUI = usejava('desktop') && feature('ShowFigureWindows');
        hWait = [];
        if useGUI
            hWait = waitbar(0, sprintf('Stage %d', stage));
        end

        % Broadcast large arrays to workers
        px   = parallel.pool.Constant(denseX);
        pcdf = parallel.pool.Constant(currentCDF);

        % Setup DataQueue for progress updates
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~) updateProgress(stage, nDense, hWait));

        % Parallel convolution integral
        parfor idx = 1:nDense
            xVal = px.Value(idx);
            if xVal > 0
                integrand = @(t) interp1(px.Value, pcdf.Value, xVal - t, 'linear', 0) .* ...
                            stage2PDF(t, aMU, vMU, aLAMBDA, vLAMBDA, w2, w2);
                newCDF(idx) = integral(integrand, 0, Inf, 'AbsTol',1e-8, 'RelTol',1e-6);
            end
            send(dq, 1);
        end

        if ~isempty(hWait) && ishandle(hWait)
            close(hWait);
        end

        % Ensure non-decreasing
        currentCDF = cummax(newCDF);
    end

    % Interpolate back to original xx
    yy = interp1(denseX, currentCDF, xx, 'linear', 'extrap');
    yy = yy + eps * randn(size(yy));
end

% Local helpers unchanged
function F = stage1CDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    F1 = cdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    F  = F1 + F2 - F1 .* F2;
    F(isnan(F)) = 0;
end

function p = stage2PDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    F1 = cdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    f1 = pdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    f2 = pdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    p  = f1 .* (1 - F2) + f2 .* (1 - F1);
    p(isnan(p)) = 0;
end

function updateProgress(stage, total, hWait)
    persistent countMap
    if isempty(countMap)
        countMap = containers.Map('KeyType','double','ValueType','double');
    end
    if ~isKey(countMap, stage)
        countMap(stage) = 0;
    end
    countMap(stage) = countMap(stage) + 1;
    frac = countMap(stage) / total;
    if ~isempty(hWait) && ishandle(hWait)
        waitbar(frac, hWait, sprintf('Stage %d: %.1f%%', stage, frac*100));
    else
        fprintf('Stage %d: %.1f%%\n', stage, frac*100);
    end
    if countMap(stage) >= total
        remove(countMap, stage);
    end
end

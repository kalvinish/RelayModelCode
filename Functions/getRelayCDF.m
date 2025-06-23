function yy = getRelayCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2)
% getRelayCDF  Relay-model CDF for one or two-stage race
%
%   yy = getRelayCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, ...)
%   computes the cumulative distribution function (CDF) of the relay model
%   with two-stage processing: stage1 then stage2.
%
%   Inputs:
%     xx       - Vector of response-time values (non-negative, increasing)
%     aMU      - Auditory IG mean
%     vMU      - Visual   IG mean
%     aLAMBDA  - Auditory IG lambda
%     vLAMBDA  - Visual   IG lambda
%     aW1      - Auditory weight for stage1
%     aW2      - Auditory weight for stage2
%     vW1      - Visual   weight for stage1
%     vW2      - Visual   weight for stage2
%
%   Output:
%     yy       - Vector of CDF values at each xx
%
    % Validate inputs
    validateattributes(xx,      {'numeric'}, {'vector','nonnegative','real','increasing'}, mfilename, 'xx', 1);
    validateattributes([aMU,vMU,aLAMBDA,vLAMBDA], {'numeric'},{'scalar','positive','real'}, mfilename);
    validateattributes([aW1,aW2,vW1,vW2], {'numeric'},{'vector','nonnegative','real','<=',1}, mfilename);

    % Preallocate output
    xx = xx(:);
    yy = zeros(size(xx));

    % Single-stage case: any weight zero => treat as one-stage race
    if any([aW1,aW2,vW1,vW2] == 0)
        yy = oneStageCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
        return;
    end

    % Two-stage relay: convolution of stage1 CDF and stage2 PDF
    % Define stage1 CDF helper
    stage1 = @(t) stage1CDF(t, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1);
    % Define stage2 PDF helper
    stage2 = @(t) stage2PDF(t, aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2);

    % Compute convolution integral for each xx
    parfor i = 1:length(xx)
        tVal = xx(i);
        if tVal > 0
            integrand = @(t) stage1(tVal - t) .* stage2(t);
            yy(i) = integral(integrand, 0, Inf, 'AbsTol',1e-6,'RelTol',1e-4);
        else
            yy(i) = 0;
        end
    end
end

%% Helper: one-stage combined CDF
function F = oneStageCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)
    % Combined race CDF: F1 + F2 - F1.*F2
    F1 = cdf('InverseGaussian', xx, aMU, aLAMBDA);
    F2 = cdf('InverseGaussian', xx, vMU, vLAMBDA);
    F  = F1 + F2 - F1 .* F2;
    F(isnan(F)) = 0;
end

%% Helper: stage1 CDF
function F = stage1CDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    % Weighted IG CDFs with race combination
    F1 = cdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    F  = F1 + F2 - F1 .* F2;
    F(isnan(F)) = 0;
end

%% Helper: stage2 PDF
function p = stage2PDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    % PDF of second-stage relay: f1*(1-F2) + f2*(1-F1)
    F1 = cdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    f1 = pdf('InverseGaussian', xx, aMU * aW, aLAMBDA * aW^2);
    f2 = pdf('InverseGaussian', xx, vMU * vW, vLAMBDA * vW^2);
    p  = f1 .* (1 - F2) + f2 .* (1 - F1);
    p(isnan(p)) = 0;
end
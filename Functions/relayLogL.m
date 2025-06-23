function [ll, pdfVals] = relayLogL(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2)
% multiLogL  Log‑likelihood for two‑stage relay model
%
%   ll = multiLogL(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2)
%
%   Inputs:
%     xx        - vector of observed response times
%     aMU, vMU  - means for A- and V-channels
%     aLAMBDA,
%     vLAMBDA   - dispersion parameters for the inverse‑Gaussian
%     aW1,vW1   - weights for stage 1
%     aW2,vW2   - weights for stage 2
%
%   Output:
%     ll        - scalar log‑likelihood: sum(log PDF(xx))

    % preallocate PDF vector
    n = numel(xx);
    pdfVals = zeros(size(xx));
    
    if (aW1==0) || (aW2==0) || (vW1==0) || (vW2==0)
        % degenerate: single-stage
        pdfVals = just1stagePDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
    else
        % true two-stage convolution
        parfor i = 1:n
            ti = xx(i);
            if ti > 0
                convFun = @(s) stage1PDF(ti - s, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1) ...
                             .* stage2PDF(s,  aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2);
                pdfVals(i) = integral(convFun, 0, Inf);
            end
        end
    end
    
    % guard against zeros
    pdfVals(pdfVals < eps) = eps;
    
    % sum log‑likelihood
    ll = sum(log(pdfVals));
end

function y = stage1PDF(x, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    % PDF of stage 1 (Inverse‑Gaussian channels A & V)
    F1 = cdf("InverseGaussian", x, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf("InverseGaussian", x, vMU * vW, vLAMBDA * vW^2);
    f1 = pdf("InverseGaussian", x, aMU * aW, aLAMBDA * aW^2);
    f2 = pdf("InverseGaussian", x, vMU * vW, vLAMBDA * vW^2);
    y  = f1 .* (1 - F2) + f2 .* (1 - F1);
    y(isnan(y)) = 0;
end

function y = stage2PDF(x, aMU, vMU, aLAMBDA, vLAMBDA, aW, vW)
    % PDF of stage 2 (same form as stage1PDF)
    F1 = cdf("InverseGaussian", x, aMU * aW, aLAMBDA * aW^2);
    F2 = cdf("InverseGaussian", x, vMU * vW, vLAMBDA * vW^2);
    f1 = pdf("InverseGaussian", x, aMU * aW, aLAMBDA * aW^2);
    f2 = pdf("InverseGaussian", x, vMU * vW, vLAMBDA * vW^2);
    y  = f1 .* (1 - F2) + f2 .* (1 - F1);
    y(isnan(y)) = 0;
end

function y = just1stagePDF(x, aMU, vMU, aLAMBDA, vLAMBDA)
    % Single-stage PDF (competing channels)
    F1 = cdf("InverseGaussian", x, aMU, aLAMBDA);
    F2 = cdf("InverseGaussian", x, vMU, vLAMBDA);
    f1 = pdf("InverseGaussian", x, aMU, aLAMBDA);
    f2 = pdf("InverseGaussian", x, vMU, vLAMBDA);
    y  = f1 .* (1 - F2) + f2 .* (1 - F1);
    y(isnan(y)) = 0;
end

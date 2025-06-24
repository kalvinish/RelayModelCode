function yy = getRaabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)
%GETRAABCDF Compute the Raab cumulative distribution function (CDF)
%   yy = GETRAABCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA) returns the combined
%   CDF according to Raab's integrated race model by combining two
%   unit CDFs while accounting for their joint probability.
%
% Inputs:
%   xx      - Numeric vector or array of time points at which to evaluate.
%   aMU     - Location parameter (mu) for the first uniCDF process.
%   vMU     - Location parameter (mu) for the second uniCDF process.
%   aLAMBDA - Scale/rate parameter (lambda) for the first uniCDF process.
%   vLAMBDA - Scale/rate parameter (lambda) for the second uniCDF process.
%
% Output:
%   yy      - Numeric vector or array of the same size as xx containing the
%             Raab CDF values, computed as F1 + F2 - F1.*F2.
%
% Example:
%   xx = 0:0.1:5;
%   yy = getRaabCDF(xx, 2.0, 3.0, 1.5, 2.0);
%   % yy contains values of F1 + F2 - F1.*F2
%
% See also uniCDF, getGriceCDF, getMillerCDF

    % Evaluate the first unit CDF at each time point
    F1 = getUniCDF(xx, aMU, aLAMBDA);

    % Evaluate the second unit CDF at each time point
    F2 = getUniCDF(xx, vMU, vLAMBDA);

    % Compute Raab CDF using inclusion-exclusion principle:
    % P(A or B) = P(A) + P(B) - P(A)P(B)
    yy = F1 + F2 - F1 .* F2;
end

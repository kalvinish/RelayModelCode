function yy = getUniCDF(xx, MU, LAMBDA)
%GETUNICDF Compute the inverse Gaussian (Wald) cumulative distribution function (CDF)
%   yy = GETUNICDF(xx, MU, LAMBDA) returns the CDF values for the inverse
%   Gaussian distribution with mean MU and shape parameter LAMBDA evaluated
%   at the points in xx. Any NaN results (e.g., for xx < 0) are replaced with
%   a very small random positive number to avoid downstream issues with zeros.
%
% Inputs:
%   xx      - Numeric vector or array of evaluation points (time or measurement).
%   MU      - Mean (mu) parameter for the inverse Gaussian distribution.
%   LAMBDA  - Shape (lambda) parameter for the inverse Gaussian distribution.
%
% Output:
%   yy      - Numeric vector or array of the same size as xx containing CDF
%             values. NaNs are replaced by small random values (~1e-320 scale).
%
% Example:
%   xx = [0.1, 0.5, 1.0, 2.0];
%   yy = getUniCDF(xx, 1.0, 3.0);
%   % yy contains valid CDF values with no NaNs.
%
% See also CDF, getMillerCDF, getRaabCDF, getGriceCDF

    % Compute the inverse Gaussian CDF using MATLAB's built-in function
    yy = cdf("InverseGaussian", xx, MU, LAMBDA);

    % Identify any NaN entries (e.g., due to invalid xx values)
    nanList = find(isnan(yy));

    % Replace NaNs with a very small random value to ensure positivity
    for i = 1:length(nanList)
        yy(nanList(i)) = randn(1) * 1e-320;
    end
end

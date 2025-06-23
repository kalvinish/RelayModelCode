function violation = getViolation_fromCDF(xx, empCDF, millerCDF)
% getViolation_fromCDF Compute Miller-bound violations from CDF data
%
%   violation = getViolation_fromCDF(xx, empCDF, millerCDF) calculates the
%   total amount by which the empirical CDF (empCDF) exceeds the theoretical
%   Miller bound CDF (millerCDF) across response times xx, using the
%   trapezoidal rule for integration.
%
%   Inputs:
%     xx         - Numeric vector of response times (monotonically increasing)
%     empCDF     - Empirical cumulative distribution values (same size as xx)
%     millerCDF  - Miller-bound CDF values (same size as xx)
%
%   Output:
%     violation  - Scalar representing the area where empCDF > millerCDF
%
%   Example:
%     xx = linspace(0,1,200);
%     emp = normcdf(xx, 0.5, 0.1);
%     miller = min(normcdf(xx, 0.4, 0.2) + normcdf(xx, 0.6, 0.2), 1);
%     v = getViolation_fromCDF(xx, emp, miller);

    % Input validation
    validateattributes(xx,        {'numeric'}, {'vector','real','increasing'}, mfilename, 'xx', 1);
    validateattributes(empCDF,    {'numeric'}, {'vector','real','numel',numel(xx)}, mfilename, 'empCDF', 2);
    validateattributes(millerCDF, {'numeric'}, {'vector','real','numel',numel(xx)}, mfilename, 'millerCDF', 3);

    % Ensure column orientation for consistency
    xx         = xx(:);
    empCDF     = empCDF(:);
    millerCDF  = millerCDF(:);

    % Calculate difference: empirical minus Miller bound
    diffCDF = empCDF - millerCDF;

    % Retain only positive violations
    diffCDF(diffCDF < 0) = 0;

    % Integrate positive differences using trapezoidal method
    violation = trapz(xx, diffCDF);
end

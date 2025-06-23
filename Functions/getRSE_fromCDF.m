function gain = getRSE_fromCDF(xx, empCDF, griceCDF)
% getRSE_fromCDF Calculate redundancy gain (RSE) from empirical and race model CDFs
%
%   gain = getRSE_fromCDF(xx, empCDF, griceCDF) computes the area where the
%   empirical cumulative distribution (empCDF) exceeds the race-model (Grice)
%   CDF (griceCDF) over the response-time vector xx using the trapezoidal rule.
%
%   Inputs:
%     xx        - Numeric vector of response times (monotonically increasing)
%     empCDF    - Empirical CDF values at each xx (same size as xx)
%     griceCDF  - Race-model (independent-channels) CDF at each xx (same size)
%
%   Output:
%     gain      - Scalar redundancy gain (area where empCDF > griceCDF)
%
%   Example:
%     xx = linspace(0,2,100);
%     emp = normcdf(xx, 1, 0.2);
%     race = normcdf(xx, 1.1, 0.3);
%     g = getRSE_fromCDF(xx, emp, race);

    % Input validation
    validateattributes(xx,       {'numeric'}, {'vector','real','increasing'}, mfilename, 'xx', 1);
    validateattributes(empCDF,   {'numeric'}, {'vector','real','numel',numel(xx)}, mfilename, 'empCDF', 2);
    validateattributes(griceCDF, {'numeric'}, {'vector','real','numel',numel(xx)}, mfilename, 'griceCDF', 3);

    % Ensure column vectors
    xx       = xx(:);
    empCDF   = empCDF(:);
    griceCDF = griceCDF(:);

    % Compute pointwise difference
    diffCDF = empCDF - griceCDF;

    % Keep only positive differences (redundancy region)
    diffCDF(diffCDF < 0) = 0;

    % Integrate using trapezoidal rule
    gain = trapz(xx, diffCDF);
end

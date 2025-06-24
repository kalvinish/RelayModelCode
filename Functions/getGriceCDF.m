function griceCDF = getGriceCDF(uni1CDF, uni2CDF)
%GETGRICECDF Compute the Grice cumulative distribution function (CDF)
%   griceCDF = GETGRICECDF(uni1CDF, uni2CDF) returns the element-wise maximum
%   of two input cumulative distribution functions uni1CDF and uni2CDF.
%   This represents the Grice CDF, selecting at each point the higher
%   cumulative probability between the two distributions.
%
% Inputs:
%   uni1CDF - Numeric vector or array representing the first CDF.
%   uni2CDF - Numeric vector or array representing the second CDF.
%             Must be the same size as uni1CDF.
%
% Output:
%   griceCDF - Numeric vector or array of the same size as the inputs,
%              containing the element-wise maximum values of uni1CDF
%              and uni2CDF.
%
% Example:
%   cdf1 = [0.1, 0.4, 0.8, 1.0];
%   cdf2 = [0.2, 0.3, 0.7, 1.0];
%   gCDF = getGriceCDF(cdf1, cdf2);
%   % gCDF = [0.2, 0.4, 0.8, 1.0]
%
% See also MAX, CDF

    % Ensure inputs are of the same size
    assert(isequal(size(uni1CDF), size(uni2CDF)), ...
        'Inputs uni1CDF and uni2CDF must be the same size.');

    % Compute the Grice CDF by taking the element-wise maximum
    griceCDF = max(uni1CDF, uni2CDF);
end

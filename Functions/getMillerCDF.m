function yy = getMillerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)
%GETMILLERCDF Compute the Miller cumulative distribution function (CDF)
%   yy = GETMILLERCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA) returns the combined
%   CDF according to Miller's race model inequality by summing two unit CDFs and
%   capping the result at 1 to maintain valid probability values.
%
% Inputs:
%   xx      - Numeric vector or array of time points to evaluate the CDF.
%   aMU     - Location parameter (mu) for the first uniCDF process.
%   vMU     - Location parameter (mu) for the second uniCDF process.
%   aLAMBDA - Scale/rate parameter (lambda) for the first uniCDF process.
%   vLAMBDA - Scale/rate parameter (lambda) for the second uniCDF process.
%
% Output:
%   yy      - Numeric vector or array of the same size as xx containing the
%             Miller CDF values, capped at a maximum of 1.
%
% Example:
%   xx = 0:0.1:5;
%   yy = getMillerCDF(xx, 2.0, 3.0, 1.5, 2.0);
%   % yy contains values of F1+F2 capped at 1
%
% See also uniCDF, getGriceCDF

    % Evaluate the first unit CDF at each time point
    F1 = getUniCDF(xx, aMU, aLAMBDA);

    % Evaluate the second unit CDF at each time point
    F2 = getUniCDF(xx, vMU, vLAMBDA);

    % Sum the two CDFs according to Miller's race model
    yy = F1 + F2;

    % Cap probabilities at 1 to ensure valid CDF values
    yy(yy > 1) = 1;
end

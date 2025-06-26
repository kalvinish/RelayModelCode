function griceCDF = getGriceCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)
%GETGRICECDF  Compute the Grice bound CDF from two unisensory CDFs
%
%   griceCDF = GETGRICECDF(xx, aMU, vMU, aLAMBDA, vLAMBDA) computes the
%   Grice cumulative distribution function (CDF) for audiovisual reaction
%   times by taking the element-wise maximum of the auditory and visual
%   unisensory CDFs at each time point in xx.
%
%   Inputs:
%     xx       — Numeric vector of timepoints at which to evaluate the CDFs
%     aMU      — Mean (μ) parameter for the auditory unisensory distribution
%     vMU      — Mean (μ) parameter for the visual  unisensory distribution
%     aLAMBDA  — Rate (λ) parameter for the auditory unisensory distribution
%     vLAMBDA  — Rate (λ) parameter for the visual  unisensory distribution
%
%   Output:
%     griceCDF — Numeric vector of the same size as xx, containing the
%                Grice bound CDF, i.e. the element-wise maximum of the
%                auditory and visual unisensory CDFs.
%
%   The function internally calls getUniCDF(xx, mu, lambda) to compute each
%   unisensory CDF:
%       uni1CDF = getUniCDF(xx, aMU, aLAMBDA);
%       uni2CDF = getUniCDF(xx, vMU, vLAMBDA);
%   and then does:
%       griceCDF = max(uni1CDF, uni2CDF);
%
%   Example:
%     xx = 0:10:100; 
%     audCDF = getUniCDF(xx, 200, 0.01);
%     visCDF = getUniCDF(xx, 250, 0.008);
%     gCDF = getGriceCDF(xx, 200, 250, 0.01, 0.008);
%
%   See also getUniCDF, getRaabCDF

    % Compute each unisensory CDF
    uniAud = getUniCDF(xx, aMU, aLAMBDA);
    uniVis = getUniCDF(xx, vMU, vLAMBDA);

    % Verify same dimensions
    assert(isequal(size(uniAud), size(uniVis)), ...
        'getGriceCDF:SizeMismatch', ...
        'Auditory and visual CDF outputs must be the same size.');

    % Take element-wise maximum to form the Grice bound
    griceCDF = max(uniAud, uniVis);
end

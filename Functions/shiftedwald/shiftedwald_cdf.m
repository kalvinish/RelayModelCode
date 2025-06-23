function F = shiftedwald_cdf(t, mu, lambda, shiftDist, shiftParams)
%MIXURESHIFTEDWALDCDF  CDF of T = X + C, where X ~ Wald(mu, lambda) and
% C ~ shiftDist with parameters shiftParams.
%
%   F = MIXURESHIFTEDWALDCDF(t, mu, lambda, shiftDist, shiftParams)
%
%   Inputs:
%       t           - scalar or vector of points at which to evaluate the CDF.
%       mu          - Wald distribution parameter (mu > 0).
%       lambda      - Wald distribution parameter (lambda > 0).
%       shiftDist   - string specifying distribution of C:
%                       'uniform', 'normal', or 'exponential'.
%       shiftParams - parameters of the C distribution, e.g.
%                       * 'uniform': [a b]
%                       * 'normal':  [mu_c, sigma_c]
%                       * 'exponential': [lambda_c]
%
%   Output:
%       F           - array of the same size as t, giving F_T(t).
%
%   Example:
%       % Suppose X ~ Wald(mu=1,lambda=2) and C ~ Uniform(0,1).
%       % Evaluate the mixture CDF at t = [0,0.5,1,1.5,...].
%       tVals = 0:0.5:3;
%       Fvals = mixtureShiftedWaldCDF(tVals, 1, 2, 'uniform', [0 1]);
%       disp([tVals(:) Fvals(:)]);

    % Convert t to a column vector for convenience
    t = t(:);
    n = numel(t);

    % Prepare output
    F = zeros(size(t));
    
    % -- Local function handle for the standard Wald CDF ----------------
    waldCDF = @(x) localWaldCDF(x, mu, lambda);
    
    % -- Define PDF of the shift distribution, plus integration limits --
    switch lower(shiftDist)
        case 'uniform'
            a = shiftParams(1);
            b = shiftParams(2);
            if a >= b
                error('For Uniform, shiftParams must be [a b] with a < b.');
            end
            pdfC = @(c) (1/(b-a)) * ones(size(c)); 
            cLower = a;
            cUpper = b;
            
        case 'normal'
            mu_c    = shiftParams(1);
            sigma_c = shiftParams(2);
            if sigma_c <= 0
                error('Normal distribution requires sigma_c > 0');
            end
            pdfC = @(c) (1/(sigma_c*sqrt(2*pi))) * ...
                        exp(-0.5*((c - mu_c)/sigma_c).^2);
            cLower = -Inf;
            cUpper = Inf;
            
        case 'exponential'
            lambda_c = shiftParams(1);
            if lambda_c <= 0
                error('Exponential distribution requires lambda_c > 0');
            end
            pdfC = @(c) lambda_c * exp(-lambda_c * c) .* (c>=0);
            cLower = 0;
            cUpper = Inf;
            
        otherwise
            error('Unsupported shiftDist: %s', shiftDist);
    end

    % -- Compute mixture CDF by numerical integration for each t(i) -----
    for i = 1:n
        ti = t(i);
        
        % Define integrand(c) = f_C(c) * F_Wald(t - c)
        integrand = @(c) pdfC(c) .* waldCDF(ti - c);
        
        % Perform numerical integration
        % We integrate over the support of C
        F(i) = integral(integrand, cLower, cUpper);
    end
end

% =======================================================================
function Fw = localWaldCDF(x, mu, lambda)
%LOCALWALDCDF  Wald/Inverse Gaussian CDF for x (vector ok).
%   x must be a real array; for x <= 0, the CDF is zero.

    Fw = zeros(size(x));
    
    % Indices where x > 0
    idxPos = (x > 0);
    xPos = x(idxPos);

    % Compute the Wald/IG CDF for x > 0
    A = sqrt(lambda ./ xPos) .* ( (xPos./mu) - 1 );
    B = sqrt(lambda ./ xPos) .* ( (xPos./mu) + 1 );

    PhiA = 0.5 * (1 + erf(A./sqrt(2)));
    PhiB = 0.5 * (1 + erf(-B./sqrt(2)));  % = Phi(-B)

    Fwald = PhiA + exp((2*lambda)./mu) .* PhiB;
    
    % Assign results
    Fw(idxPos) = Fwald;
    
    % For x <= 0, remains 0
end

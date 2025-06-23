% ExGaussian Unisensory (just a standard exgaussian distribution)
function yy = exgauss_cdf(xx, mu, sigma, tau)
    % Vectorized ex-Gaussian CDF
    % xx, mu, sigma, tau can be scalars or arrays of matching size
    
    % Standardize
    z1 = (xx - mu)./sigma;
    z2 = z1 - sigma./tau;
    
    % First normal CDF
    term1 = normcdf(z1);
    
    % Exponential factor times second normal CDF
    exp_factor = exp((sigma.^2)/(2*tau^2) - (xx - mu)./tau);
    term2 = exp_factor .* normcdf(z2);
    
    % Final ex-Gaussian CDF
    yy = term1 - term2;
end
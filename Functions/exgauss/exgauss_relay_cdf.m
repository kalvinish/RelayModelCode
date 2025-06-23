% Relay between exgaussian distributions (CDF)
function yy = exgauss_relay_cdf(xx, mu, sigma, tau)
    parfor ii=1:length(xx)
        if xx(ii)>0
            fun = @(t) exponential_race_cdf(xx(ii)-t, tau) .*...
                       gaussian_race_pdf(t, mu, sigma);

            yy(ii) = integral( fun, -Inf, Inf );
        end
    end 
end

% Race between two exponential distributions (CDF)
function yy = exponential_race_cdf(xx, tau)
    F1 = expcdf(xx, tau(1));
    F2 = expcdf(xx, tau(2));

    yy = F1 + F2 - (F1 .* F2);
end

% Race between two Gaussian distributions (PDF)
function yy = gaussian_race_pdf(xx, mu, sigma)
    F1 = normcdf(xx, mu(1), sigma(1));
    F2 = normcdf(xx, mu(2), sigma(2));
    f1 = normpdf(xx, mu(1), sigma(1));
    f2 = normpdf(xx, mu(2), sigma(2));
    
    yy = f1 .* (1 - F2) + f2 .* (1 - F1);
    
    yy(isnan(yy)) = 0;
end
% EcGaussian Race
function yy = exgauss_race_cdf(xx, mu, sigma, tau)
    F1 = exgauss_cdf(xx, mu(1), sigma(1), tau(1));
    F2 = exgauss_cdf(xx, mu(2), sigma(2), tau(2));

    yy = F1 + F2 - (F1 .* F2);
end
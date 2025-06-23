function [mu_opt, lambda_opt, rmse, ks_stat, h, p, boot_pvalue] = fit_inverse_gaussian_cdf(data)
    % Fits an Inverse Gaussian CDF by directly minimising RMSE.
    %
    % Inputs:
    %   data - N×2: [x, emp_CDF(x)]
    %
    % Outputs:
    %   mu_opt, lambda_opt — fitted parameters
    %   rmse      — root-mean-square error of the fit
    %   ks_stat,h,p — KS test results vs. fitted CDF
    %   boot_pvalue — bootstrap KS p-value

    % Extract and dedupe
    x           = data(:,1);
    F_empirical = data(:,2);
    [x, idx]    = unique(x);
    F_empirical = F_empirical(idx);
    n           = numel(x);

    % Initial guess
    mu0     = mean(x);
    lambda0 = 1;
    params0 = [mu0, lambda0];

    % Model CDF handle
    modelCDF = @(params) cdf('InverseGaussian', x, params(1), params(2));

    % Objective: RMSE
    objfun = @(params) sqrt( mean( (F_empirical - modelCDF(params)).^2 ) );

    % Constrain mu>0, lambda>0
    lb = [eps, eps];
    ub = [];

    opts = optimoptions('fmincon', ...
               'Display','off', ...
               'TolFun',1e-10, ...
               'TolX',1e-10);

    % Run the minimisation
    [params_opt, rmse] = fmincon(objfun, params0, [],[],[],[], lb, ub, [], opts);

    mu_opt     = params_opt(1);
    lambda_opt = params_opt(2);

    % Compute fitted CDF and do KS test
    F_fit       = modelCDF(params_opt);
    [h, p, ks_stat] = kstest(x, 'CDF', [x, F_fit]);

    % Bootstrap for KS p-value
    B             = 10000;
    ks_bootstrap  = zeros(B,1);
    for b = 1:B
        u       = rand(n,1);
        x_sim   = icdf('InverseGaussian', u, mu_opt, lambda_opt);
        x_sim   = min(max(x_sim, min(x)), max(x));
        [~, ~, ks_bootstrap(b)] = kstest(x_sim, 'CDF', [x, F_fit]);
    end
    boot_pvalue = mean(ks_bootstrap >= ks_stat);
end
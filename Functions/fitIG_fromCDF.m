function [mu_opt, lambda_opt, rmse] = fitIG_fromCDF(data)
    % Fits an Inverse Gaussian CDF by directly minimising RMSE.
    %
    % Inputs:
    %   data - N×2: [x, emp_CDF(x)]
    %
    % Outputs:
    %   mu_opt, lambda_opt — fitted parameters
    %   rmse      — root-mean-square error of the fit

    % Extract and dedupe
    x           = data(:,1);
    F_empirical = data(:,2);
    [x, idx]    = unique(x);
    F_empirical = F_empirical(idx);

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

end
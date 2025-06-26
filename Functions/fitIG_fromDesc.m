function [mu, lambda, rmse] = fitIG_fromDesc(target_mean, target_se, n)
    % Fit an Inverse Gaussian so as to minimise the RMSE between
    % (mean, SD, median) and their targets.
    %
    % Inputs:
    %   target_mean   – desired mean
    %   target_se     – desired standard error
    %   n             – sample size used to compute SE
    %
    % Outputs:
    %   mu, lambda – fitted parameters
    %   rmse        – the achieved root‐mean‐square‐error

    % Convert SE → SD
    target_sd = target_se * sqrt(n);

    % Initial guesses (method‐of‐moments for mu,lambda)
    mu0     = target_mean;
    lambda0 = target_mean^3 / target_sd^2;
    x0      = [mu0, lambda0];

    % Objective: RMSE of the three summary stats
    objfun = @(params) sqrt( mean( residuals(params, target_mean, target_sd).^2 ) );

    % Constrain mu>0, lambda>0
    lb = [eps, eps];
    ub = [];

    opts = optimoptions('fmincon', ...
                       'Display','off', ...
                       'TolX',1e-10, ...
                       'TolFun',1e-10);

    % Run the minimiser
    params_opt = fmincon(objfun, x0, [], [], [], [], lb, ub, [], opts);

    % Unpack and compute final RMSE
    mu     = params_opt(1);
    lambda = params_opt(2);
    rmse   = objfun(params_opt);
end

function res = residuals(params, tgt_m, tgt_s)
    % Compute the three raw errors
    mu     = params(1);
    lambda = params(2);

    % Enforce positivity
    if mu<=0 || lambda<=0
        res = [Inf;Inf;Inf];
        return;
    end

    % 1) mean error
    err_mean = mu - tgt_m;

    % 2) SD error from var = mu^3/lambda
    err_sd = sqrt(mu^3/lambda) - tgt_s;

    % res = [err_mean; err_sd; err_med];
    res = [err_mean; err_sd];
end

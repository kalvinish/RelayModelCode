function yy = uniRND(n, MU, LAMBDA, sortFlag)
    % uniRND generates n random numbers from the Inverse Gaussian distribution
    % with parameters MU and LAMBDA.
    % If sortFlag is provided and set to 1, the output is sorted in ascending order.
    % By default, sortFlag is 0, and the output is unsorted.

    % Check if sortFlag is provided; if not, default to 0 (unsorted)
    if nargin < 4
        sortFlag = 0;
    end

    % Generate random samples from the Inverse Gaussian distribution
    yy = random("InverseGaussian", MU, LAMBDA, n, 1);

    % Sort the output if sortFlag is 1
    if sortFlag == 1
        yy = sort(yy, 'ascend');
    end
end
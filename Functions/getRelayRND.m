function rt = getRelayRND(n, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2, sortFlag)
    % getRelayRND generates reaction times based on auditory and visual parameters.
    % If sortFlag is provided and set to 1, the output rt is sorted in ascending order.
    % By default, sortFlag is 0, and the output rt is unsorted.
    
    % Check if sortFlag is provided; if not, default to 0 (unsorted)
    if nargin < 10  % There are 9 input arguments before sortFlag
        sortFlag = 0;
    end

    % Check if any of the weights are zero
    if (aW1 == 0) || (aW2 == 0) || (vW1 == 0) || (vW2 == 0)
        % Generate random samples when any weight is zero
        rt = min([random('inversegaussian', aMU * (aW1 + aW2), aLAMBDA * ((aW1 + aW2)^2), n, 1), ...
                  random('inversegaussian', vMU * (vW1 + vW2), vLAMBDA * ((vW1 + vW2)^2), n, 1)], [], 2);
    else
        % Generate random samples for the first component
        first = min([random('inversegaussian', aMU * aW1, aLAMBDA * aW1^2, n, 1), ...
                     random('inversegaussian', vMU * vW1, vLAMBDA * vW1^2, n, 1)], [], 2);
        % Generate random samples for the second component
        second = min([random('inversegaussian', aMU * aW2, aLAMBDA * aW2^2, n, 1), ...
                      random('inversegaussian', vMU * vW2, vLAMBDA * vW2^2, n, 1)], [], 2);
        % Sum the two components to get the reaction time
        rt = first + second;
    end

    % Sort the output if sortFlag is 1
    if sortFlag == 1
        rt = sort(rt, 'ascend');
    end
end
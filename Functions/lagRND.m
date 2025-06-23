% function rt = lagRND( n, mu1, mu2, lambda1, lambda2, w1_1, w1_2, w2_1, w2_2, lag, sortFlag )
%     % lagRND generates reaction times based on the specified parameters.
%     % If lag is positive, it is added to dist1 in both races.
%     % If lag is negative, it is added to dist2 in both races.
%     % Any remaining lag after race 1 is added to the same distribution in race 2.
%     % If sortFlag is provided and set to 1, the output rt is sorted in ascending order.
%     % By default, sortFlag is 0, and the output rt is unsorted.
% 
%     % Check if sortFlag is provided; if not, default to 0 (unsorted)
%     if nargin < 11  % There are 10 input arguments before sortFlag
%         sortFlag = 0;
%     end
% 
%     if (w1_1 > 0) && (w2_1 > 0)
%         % Generate random samples from the Inverse Gaussian distributions
%         r_dist1_1 = random( 'InverseGaussian', mu1 * w1_1, lambda1 * w1_1^2, [n, 1] );
%         r_dist1_2 = random( 'InverseGaussian', mu1 * w1_2, lambda1 * w1_2^2, [n, 1] );
%         r_dist2_1 = random( 'InverseGaussian', mu2 * w2_1, lambda2 * w2_1^2, [n, 1] );
%         r_dist2_2 = random( 'InverseGaussian', mu2 * w2_2, lambda2 * w2_2^2, [n, 1] );
% 
%         if lag >= 0
%             % Positive lag: lag is added to dist1 in both races
%             adjusted_r_dist1_1 = r_dist1_1 + lag;
%             adjusted_r_dist2_1 = r_dist2_1;
% 
%             % Compute the winner of race 1
%             r1win = min( adjusted_r_dist1_1, adjusted_r_dist2_1 );
% 
%             % Remaining lag
%             lagRemaining = max( lag - r1win, 0 );
% 
%             % Adjust dist1 in race 2 with any remaining lag
%             adjusted_r_dist1_2 = r_dist1_2 + lagRemaining;
%             adjusted_r_dist2_2 = r_dist2_2;
%         else
%             % Negative lag: lag is added to dist2 in both races
%             adjusted_r_dist1_1 = r_dist1_1;
%             adjusted_r_dist2_1 = r_dist2_1 + abs( lag );  % lag is negative
% 
%             % Compute the winner of race 1
%             r1win = min( adjusted_r_dist1_1, adjusted_r_dist2_1 );
% 
%             % Remaining lag
%             lagRemaining = max( abs( lag ) - r1win, 0);
% 
%             % Adjust dist2 in race 2 with any remaining lag
%             adjusted_r_dist1_2 = r_dist1_2;
%             adjusted_r_dist2_2 = r_dist2_2 + lagRemaining;  % lagRemaining is negative or zero
%         end
% 
%         % Compute the winner of race 2
%         r2win = min( adjusted_r_dist1_2, adjusted_r_dist2_2 );
% 
%         % Compute the full reaction time
%         rt = r1win + r2win;
% 
%     else
% 
%         if lag < 0
%             lag1 = 0;
%             lag2 = abs( lag );
%         else
%             lag1 = abs( lag );
%             lag2 = 0;
%         end
% 
%         r_dist1 = random( 'InverseGaussian', mu1, lambda1, [n, 1] ) + lag1;
%         r_dist2 = random( 'InverseGaussian', mu2, lambda2, [n, 1] ) + lag2;
% 
%         rt = min( r_dist1, r_dist2 );
% 
%     end
% 
%     % Sort the output if sortFlag is 1
%     if sortFlag == 1
%         rt = sort( rt, 'ascend' );
%     end
% 
% end

function rt = lagRND(n, mu1, mu2, lambda1, lambda2, ...
                     w1_1, w1_2, w2_1, w2_2, lag, sortFlag)
% Fast vectorized relay-model RT generator
%
% rt = lagRND(n,mu1,mu2,λ1,λ2,w1_1,w1_2,w2_1,w2_2,lag,sortFlag)
%  - n          : # of samples
%  - mu*,lambda*: IG params
%  - w*         : mixture weights
%  - lag        : scalar (can be negative)
%  - sortFlag   : (opt) 1→sort ascending
%
% Vectorized version avoids any per-sample branching.

    if nargin<11
        sortFlag = false;
    end

    % 1) Precompute lag offsets
    lagPos = max(lag,0);
    lagNeg = max(-lag,0);

    % Shortcut: if either race uses zero weight, fall back to simple min
    if w1_1<=0 || w2_1<=0
        % unimodal races
        a = random('InverseGaussian', mu1, lambda1, [n,1]) + lagPos;
        b = random('InverseGaussian', mu2, lambda2, [n,1]) + lagNeg;
        rt = min(a,b);
    else
        % 2) One-time draws for both races
        r11 = random('InverseGaussian', mu1*w1_1, lambda1*w1_1^2, [n,1]);
        r12 = random('InverseGaussian', mu1*w1_2, lambda1*w1_2^2, [n,1]);
        r21 = random('InverseGaussian', mu2*w2_1, lambda2*w2_1^2, [n,1]);
        r22 = random('InverseGaussian', mu2*w2_2, lambda2*w2_2^2, [n,1]);

        % 3) First race with lag
        a1     = r11 + lagPos;
        b1     = r21 + lagNeg;
        r1win  = min(a1, b1);
        % remaining lag is always abs(lag) minus win-time
        lagRem = max(max(lagPos,lagNeg) - r1win, 0);

        % 4) Second race: apply remaining lag on same channel
        %    +ve lag→channel1, -ve lag→channel2
        isPos  = lag>=0;
        a2     = r12 + (lagRem * isPos);
        b2     = r22 + (lagRem * ~isPos);

        % 5) Sum winners
        r2win = min(a2, b2);
        rt    = r1win + r2win;
    end

    % 6) Optional sorting
    if sortFlag
        rt = sort(rt);
    end
end

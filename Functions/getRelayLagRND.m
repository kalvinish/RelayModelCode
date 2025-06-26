function rt = getRelayLagRND(n, mu1, mu2, lambda1, lambda2, ...
                     w1_1, w1_2, w2_1, w2_2, lag, sortFlag)
% Fast vectorized relay-model RT generator
%
% rt = getRelayLagRND(n,mu1,mu2,λ1,λ2,w1_1,w1_2,w2_1,w2_2,lag,sortFlag)
%   - lag may now be -Inf or +Inf without producing NaNs.

    if nargin<11
        sortFlag = false;
    end

    % 1) Precompute lag offsets
    lagPos = max(lag,0);
    lagNeg = max(-lag,0);

    %  shortcut for purely unimodal races
    if w1_1<=0 || w2_1<=0
        a = random('InverseGaussian',mu1,lambda1,[n,1]) + lagPos;
        b = random('InverseGaussian',mu2,lambda2,[n,1]) + lagNeg;
        rt = min(a,b);
        if sortFlag, rt = sort(rt); end
        return
    end

    % 2) Draw all relay samples
    r11 = random('InverseGaussian', mu1*w1_1, lambda1*w1_1^2, [n,1]);
    r12 = random('InverseGaussian', mu1*w1_2, lambda1*w1_2^2, [n,1]);
    r21 = random('InverseGaussian', mu2*w2_1, lambda2*w2_1^2, [n,1]);
    r22 = random('InverseGaussian', mu2*w2_2, lambda2*w2_2^2, [n,1]);

    % 3) First race with static lag
    a1    = r11 + lagPos;
    b1    = r21 + lagNeg;
    r1win = min(a1, b1);

    % 4) Remaining lag after winner
    %    if lag was infinite, lagRem will be Inf for all draws
    lagRem = max(max(lagPos,lagNeg) - r1win, 0);

    % 5) Second race: add remaining lag only to the same channel that was
    %    lagged initially
    if lag >= 0
        % positive lag: channel 1 was initially delayed
        a2 = r12 + lagRem;
        b2 = r22;
    else
        % negative lag: channel 2 was initially delayed
        a2 = r12;
        b2 = r22 + lagRem;
    end

    % 6) Sum the two winning times
    r2win = min(a2, b2);
    rt    = r1win + r2win;

    % 7) Optional sorting
    if sortFlag
        rt = sort(rt);
    end
end

function rt = getRelayLagRND(n, mu1, mu2, lambda1, lambda2, ...
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

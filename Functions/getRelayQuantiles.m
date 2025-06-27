function qs = getRelayQuantiles(ps, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2)
  % ps: vector of target CDF values
  qs = zeros(size(ps));
  
  % Define an anonymous function that, for a fixed p, gives the squared error.
  function err = cdfErr(x, p)
    val = getRelayCDF(x, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2);
    err = (val - p).^2;
  end

  % Loop over each p
  for i = 1:numel(ps)
    p = ps(i);
    % initial guess
    x0 = aMU;  
    % call fminsearch with a small tolerance
    opts = optimset('TolX',1e-6,'TolFun',1e-8);
    qs(i) = fminsearch(@(x) cdfErr(x,p), x0, opts);
  end
end


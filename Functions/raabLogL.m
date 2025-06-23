function [ll, pdfVals] = raabLogL(xx, aMU, vMU, aLAMBDA, vLAMBDA)

    % Single-stage PDF (competing channels)
    F1 = cdf("InverseGaussian", xx, aMU, aLAMBDA);
    F2 = cdf("InverseGaussian", xx, vMU, vLAMBDA);
    f1 = pdf("InverseGaussian", xx, aMU, aLAMBDA);
    f2 = pdf("InverseGaussian", xx, vMU, vLAMBDA);

    pdfVals  = f1 .* (1 - F2) + f2 .* (1 - F1);
    pdfVals(isnan(pdfVals)) = 0;

    % guard against zeros
    pdfVals(pdfVals < eps) = eps;
    
    % sum log‑likelihood
    ll = sum(log(pdfVals));

end
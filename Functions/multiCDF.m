function yy = multiCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, aW2, vW1, vW2)

yy = zeros( size(xx) );

if (aW1 == 0) || (aW2 == 0) || (vW1 == 0) || (vW2 == 0) 
    yy = just1stageCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA);
else
    parfor ii=1:length(xx)
        if xx(ii)>0
            fun = @(t) stage1CDF(xx(ii)-t, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1) .* ...
                       stage2PDF(t, aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2);
    
            yy(ii) = integral( fun, 0, Inf );
        end
    end 
end

end

%%
function yy = stage1CDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW1, vW1)

F1 = cdf("InverseGaussian", xx, aMU*aW1, aLAMBDA*aW1^2);
F2 = cdf("InverseGaussian", xx, vMU*vW1, vLAMBDA*vW1^2);

yy = F1 + F2 - (F1 .* F2);

yy(isnan(yy)) = 0;

end

%%
function yy = stage2PDF(xx, aMU, vMU, aLAMBDA, vLAMBDA, aW2, vW2)

F1 = cdf("InverseGaussian", xx, aMU*aW2, aLAMBDA*aW2^2);
F2 = cdf("InverseGaussian", xx, vMU*vW2, vLAMBDA*vW2^2);
f1 = pdf("InverseGaussian", xx, aMU*aW2, aLAMBDA*aW2^2);
f2 = pdf("InverseGaussian", xx, vMU*vW2, vLAMBDA*vW2^2);

yy = f1 .* (1 - F2) + f2 .* (1 - F1);

yy(isnan(yy)) = 0;

end

%%
function yy = just1stageCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)

F1 = cdf("InverseGaussian", xx, aMU, aLAMBDA);
F2 = cdf("InverseGaussian", xx, vMU, vLAMBDA);

yy = F1 + F2 - (F1 .* F2);

yy(isnan(yy)) = 0;

end
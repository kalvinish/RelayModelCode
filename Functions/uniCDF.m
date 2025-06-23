function yy = uniCDF(xx, MU, LAMBDA)

yy = cdf("InverseGaussian", xx, MU, LAMBDA);

nanList = find(isnan(yy));

for i = 1:length(nanList)
    yy(nanList(i)) = randn(1)*1e-320;
end

end
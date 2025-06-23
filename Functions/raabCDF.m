function yy = raabCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)

F1 = uniCDF(xx, aMU, aLAMBDA);
F2 = uniCDF(xx, vMU, vLAMBDA);

yy = F1 + F2 - (F1 .* F2);

yy(isnan(yy)) = 0;

end
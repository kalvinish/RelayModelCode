function yy = getMillerCDF(xx, aMU, vMU, aLAMBDA, vLAMBDA)

F1 = uniCDF(xx, aMU, aLAMBDA);
F2 = uniCDF(xx, vMU, vLAMBDA);

yy = F1 + F2;

yy(yy>1) = 1;

end
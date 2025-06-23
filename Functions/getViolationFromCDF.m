function violation = getViolationFromCDF(xx, empCDF, millerCDF)

miller_difference = empCDF - millerCDF;
miller_positiveDifference = max(0, miller_difference);

violation = trapz(xx, miller_positiveDifference);

end
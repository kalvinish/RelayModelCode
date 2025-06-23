function gain = getGainFromCDF(xx, empCDF, griceCDF)

gain_difference = empCDF - griceCDF;
gain_positiveDifference = max(0, gain_difference);

gain = trapz(xx, gain_positiveDifference);

end
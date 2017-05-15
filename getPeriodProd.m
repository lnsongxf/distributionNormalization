function indProd = getPeriodProd(numBins,workersPerBin,firmsPerBin,RP,wDist1,fDist1)
  indProd = zeros(numBins,numBins);
  for i1 = 1:numBins
    ind1 = (i1-1)*workersPerBin + 1 : i1*workersPerBin;
    for i2 = 1:numBins
      ind2 = (i2-1)*firmsPerBin + 1 : i2*firmsPerBin;
      ind2 = vec(repmat(ind2,workersPerBin/firmsPerBin,1))';
      indProd(i1,i2)        = mean(arrayfun(RP,wDist1(ind1),fDist1(ind2)));
    end
  end
end
function taken = sampleCDF(adoptionProb1,numberWanted)
  inpdf = adoptionProb1;
  cdf   = cumsum(inpdf)./sum(inpdf);
  count = 0;
  taken = false(size(adoptionProb1));
  while count < numberWanted
    bla = find(rand() < cdf,1,'first');
    if taken(bla) == false
      taken(bla) = true;
      inpdf(bla) = 0;
      cdf   = cumsum(inpdf)./sum(inpdf);
      count = count + 1;
    end
  end
end
function [wTrueRank,wEstRank,wTrueBin,wEstBin,fTrueRank,fEstRank,fTrueBin,fEstBin] = ...
    getNLSInputs(indProd,wAdoptionNum,addn,fAdoptionNum,numBins,iSet)
  [~,~,~,~,~,RD,SimO] = KL('customProd',indProd,wAdoptionNum,addn,iSet);
  wTrueRank    = RD.I.iNRRankAgg(:,1);
  wEstRank     = RD.I.iNRRankAgg(:,2);
  wTrueBin     = SimO.iNameX(wTrueRank);
  wEstBin      = RD.I.iBin;
  fTrueRank    = vec(1:fAdoptionNum);
  fEstRank     = RD.J.NROmega(:,2);
  fTrueBin     = ceil(vec(1:fAdoptionNum)/(fAdoptionNum/numBins));
  fEstBin      = RD.J.jBin;
end

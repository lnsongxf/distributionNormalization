for iDist2 = 1:3
  dCUse = distCenter2{iDist2};
  %Workers and firms involved in Period 2
  %Same setup as period 1.
  wName2   = sort(datasample(vec(1:nx),wAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,nx),dCUse(1),dCUse(2))));
  fName2   = sort(datasample(vec(1:ny),fAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,ny),dCUse(1),dCUse(2))));
  wDist2   = workerDist(wName2);
  fDist2   = firmDist(fName2);
  for iProd2 = 1:3
    RPUse = RP{iProd2};
    indProd2 = getPeriodProd(numBins,workersPerBin,firmsPerBin,RPUse,wDist2,fDist2);
    for iSet2 = 1:8
      addn     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
      if exist(['.',filesep,'Output',filesep,addn,'.mat'],'file') == 0
        try
          [wTrueRank2,wEstRank2,wTrueBin2,wEstBin2,fTrueRank2,fEstRank2,fTrueBin2,fEstBin2] = getNLSInputs(indProd2,wAdoptionNum,addn,fAdoptionNum,numBins,iSet2);
          save(['.',filesep,'Output',filesep,addn,'.mat'],'fName2','wName2','wTrueRank2','wEstRank2','wTrueBin2','wEstBin2','fTrueRank2','fEstRank2','fTrueBin2','fEstBin2','wDist2','fDist2')
        end
      end
    end
  end
end

for iProd = 1:3
  %Consider a true production function
  RPUse = RP{iProd};
  %Workers and firms involved in Period 1
  wName1   = sort(datasample(vec(1:nx),wAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,nx),distCenter1(1),distCenter1(2))));
  fName1   = sort(datasample(vec(1:ny),fAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,ny),distCenter1(1),distCenter1(2))));
  wDist1   = workerDist(wName1);
  fDist1   = firmDist(fName1);
  %Calculate the production function at the bin level
  indProd1 = getPeriodProd(numBins,workersPerBin,firmsPerBin,RPUse,wDist1,fDist1);
  for iSet = 1:8
    %File identifier
    addn     = ['Period1','_Prod',num2str(iProd),'_Dist1','_Set',num2str(iSet)];
    if exist(['.',filesep,'Output',filesep,addn,'.mat'],'file') == 0
      %Obtain the numbers of interest.
      [wTrueRank1,wEstRank1,wTrueBin1,wEstBin1,fTrueRank1,fEstRank1,fTrueBin1,fEstBin1] = getNLSInputs(indProd1,wAdoptionNum,addn,fAdoptionNum,numBins,iSet);
      %Save for later
      save(['.',filesep,'Output',filesep,addn,'.mat'],'fName1','wName1','wTrueRank1','wEstRank1','wTrueBin1','wEstBin1','fTrueRank1','fEstRank1','fTrueBin1','fEstBin1','wDist1','fDist1')
    end
  end
end
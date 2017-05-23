%Construct the mapping between distributions
for iProd = 1:3
  for iSet = 1:8
    addn     = ['Period1','_Prod',num2str(iProd),'_Dist1','_Set',num2str(iSet)];
    if exist(['.',filesep,'Output',filesep,addn,'.mat'],'file') > 0
      load(['.',filesep,'Output',filesep,addn,'.mat'])
      for iDist2 = 1:3
        dCUse = distCenter2{iDist2};
        for iProd2 = 1:3
          for iSet2 = 1:8
            addn     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
            if exist(['.',filesep,'Output',filesep,addn,'.mat'],'file') > 0
              load(['.',filesep,'Output',filesep,addn,'.mat'])
              %See how well we do on the CDF for now
              mapping{iDist2,iProd,iProd2,iSet,iSet2} = nlsFunc(fEstBin1,fName2,wEstBin1,wName2,...
                fEstBin2,fTrueBin1,wEstBin2,wTrueBin1,...
                fEstRank1,fTrueBin2,wEstRank1,wTrueBin2,...
                fEstRank2,fTrueRank1,wEstRank2,wTrueRank1,...
                fName1,fTrueRank2,wName1,wTrueRank2,...
                wAdoptionNum,fAdoptionNum,numBins,dCUse(1),dCUse(2));
            end
          end
        end
      end
    end
  end
end
save mapping.mat





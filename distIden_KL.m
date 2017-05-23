clear
addpath Source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx              = 100000;
%Each firm hires 100 workers at full capacity
firmSize        = 100;
%Which means firm population is 500
ny              = nx/firmSize;

%This is the number of workers in each period
wAdoptionNum    = 60000;
%Likewise the number of firms.
fAdoptionNum    = wAdoptionNum/firmSize;
%For calculating the production function, we need to use bins.
numBins         = 100;
workersPerBin   = wAdoptionNum/numBins;
firmsPerBin     = fAdoptionNum/numBins;

%This is the true distribution at the i j level.
workerDist = sort(rand(nx,1));
firmDist   = sort(rand(ny,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTUAL production functions. This is the object we want to recover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Consider RP (Real Production) defined on RX and RY
RP{1} = @(x,y) 0.6 + 0.4* (x^0.5 + y^0.5)^(2);
RP{2} = @(x,y) (x^2 + 2*y^2)^(1/2);
RP{3} = @(x,y) (0.4 + (x-0.4+1)*y).*double(x<=0.4) + (0.4+((x-0.4)^2+y^2)^(1/2)).*double(x>0.4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Period 1 distribution normalized to a uniform
distCenter1     = [1,1];
%Period 2 distributions can be different from Period 1
distCenter2{1}     = [1,1];
distCenter2{2}     = [2,4];
distCenter2{3}     = [4,2];

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
        [wTrueRank2,wEstRank2,wTrueBin2,wEstBin2,fTrueRank2,fEstRank2,fTrueBin2,fEstBin2] = getNLSInputs(indProd2,wAdoptionNum,addn,fAdoptionNum,numBins,iSet2);
        save(['.',filesep,'Output',filesep,addn,'.mat'],'fName2','wName2','wTrueRank2','wEstRank2','wTrueBin2','wEstBin2','fTrueRank2','fEstRank2','fTrueBin2','fEstBin2','wDist2','fDist2')
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

%Construct the mapping between distributions
for iProd = 1:3
  for iSet = 1:8
    addn     = ['Period1','_Prod',num2str(iProd),'_Dist1','_Set',num2str(iSet)];
    load(['.',filesep,'Output',filesep,addn,'.mat'])
    for iDist2 = 1:3
      dCUse = distCenter2{iDist2};
      for iProd2 = 1:3
        for iSet2 = 1:8
          addn     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
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

%Experimentation on changing parameters and distribs.
%Note that the effects of changing one single entity on variance is already
%computed since we are calculating every combination of the model. The
%question now is, can we obtain the same measures of increases in variance
%from 

%Consider all combinations of parameters to and parameters from.
for iProd = 1:3
  for iSet = 1:8
    addn1     = ['Period1','_Prod',num2str(iProd),'_Dist1_Set',num2str(iSet)];
    bla = dir(['Output\*',addn1,'.mat']);
    load(['Output\',bla(1).name]);
    P1Prod      = RD.Y.Prod;
    for iDist2 = 1:3
      for iProd2 = 1:3
        for iSet2 = 1:8
          addn2     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
          bla = dir(['Output\*',addn2,'.mat']);
          load(['Output\',bla(1).name]);
          P2Prod  = RD.Y.Prod;
          
          %What if only distribution had changed.
          %Take the production function estimated from the mapped distrib,
          %recompute using that.
          paramW    = mapping{iDist2,iProd,iProd2,iSet,iSet2}.work_esti_bins.params;
          paramF    = mapping{iDist2,iProd,iProd2,iSet,iSet2}.firm_esti_bins.params;
          EstWProd  = kumaraswamyiCDF(linspace(0,1,numBins),paramW(1),paramW(2));
          EstFProd  = kumaraswamyiCDF(linspace(0,1,numBins),paramF(1),paramF(2));

          %We need this interpolated on the productivities of period one.
          [meshF,meshW] = meshgrid(EstFProd,EstWProd);
          
          %First period distribution is normalized 
          TrueWProd     = grpstats(wDist1,SimO.iNameX);
          TrueFProd     = grpstats(fDist1,SimO.jNameY);
          [qF,qW]       = meshgrid(TrueFProd,TrueWProd);
          
          %Use the counterfactual production function.
          indProd1      = interp2(meshF,meshW,P2Prod,qF,qW);
          
          %Compute what happens if you use the second period production 
          bla = ['P2Prod_t1','P',num2str(iProd),'S',num2str(iSet),'t2',num2str(iProd2),'D',num2str(iDist2),'S',num2str(iSet)];
          getNLSInputs(indProd1,wAdoptionNum,bla,fAdoptionNum,numBins,iSet);
          
          bla = ['P2Dist_t1','P',num2str(iProd),'S',num2str(iSet),'t2',num2str(iProd2),'D',num2str(iDist2),'S',num2str(iSet)];
          %What if use second period distribution?
          indProd2 = interp2(qF,qW,P1Prod,meshF,meshW);
          getNLSInputs(indProd2,wAdoptionNum,bla,fAdoptionNum,numBins,iSet);
          
        end
      end
    end
  end
end


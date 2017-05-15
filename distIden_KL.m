clear
addpath Source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx              = 50000;
%Each firm hires 100 workers at full capacity
firmSize        = 100;
%Which means firm population is 500
ny              = nx/firmSize;

%This is the number of workers in each period
wAdoptionNum    = 30000;
%Likewise the number of firms.
fAdoptionNum    = wAdoptionNum/firmSize;
%For calculating the production function, we need to use bins.
numBins         = 50;
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
  %File identifier
  addn     = ['Prod1_',num2str(iProd),'Dist1'];
  %Obtain the numbers of interest.
  [wTrueRank1,wEstRank1,wTrueBin1,wEstBin1,fTrueRank1,fEstRank1,fTrueBin1,fEstBin1] = getNLSInputs(indProd1,wAdoptionNum,addn,fAdoptionNum,numBins);
  %Save for later
  save(['.',filesep,'Output',filesep,addn,'.mat'],'wName1','wTrueRank1','wEstRank1','wTrueBin1','wEstBin1','fTrueRank1','fEstRank1','fTrueBin1','fEstBin1')
  for iDist = 1:3
    dCUse = distCenter2{iDist};
    %Workers and firms involved in Period 2
    %Same setup as period 1.
    wName2   = sort(datasample(vec(1:nx),wAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,nx),dCUse(1),dCUse(2))));
    fName2   = sort(datasample(vec(1:ny),fAdoptionNum,'Replace',false,'Weights',betapdf(linspace(0,1,ny),dCUse(1),dCUse(2))));
    wDist2   = workerDist(wName2);
    fDist2   = firmDist(fName2);
    indProd2 = getPeriodProd(numBins,workersPerBin,firmsPerBin,RPUse,wDist2,fDist2);
    addn     = ['Prod2_',num2str(iProd),'Dist',num2str(iDist)];
    [wTrueRank2,wEstRank2,wTrueBin2,wEstBin2,fTrueRank2,fEstRank2,fTrueBin2,fEstBin2] = getNLSInputs(indProd2,wAdoptionNum,addn,fAdoptionNum,numBins);
    save(['.',filesep,'Output',filesep,addn,'.mat'],'wName2','wTrueRank2','wEstRank2','wTrueBin2','wEstBin2','fTrueRank2','fEstRank2','fTrueBin2','fEstBin2')
    %See how well we do on the CDF for now
    results{iProd,iDist} = nlsFunc(fEstBin1,fName2,wEstBin1,wName2,...
      fEstBin2,fTrueBin1,wEstBin2,wTrueBin1,...
      fEstRank1,fTrueBin2,wEstRank1,wTrueBin2,...
      fEstRank2,fTrueRank1,wEstRank2,wTrueRank1,...
      fName1,fTrueRank2,wName1,wTrueRank2,...
      wAdoptionNum,fAdoptionNum,numBins,dCUse(1),dCUse(2));
  end
end




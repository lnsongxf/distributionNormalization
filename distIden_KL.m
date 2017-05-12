clear
addpath Source

%Labels for worker and firm original names.


%Consider RP (Real Production) defined on RX and RY 
%This is the TRUE output which doesn't change over time.
rAlpha = 0.3;
RP = @(rx,ry)  0.5 + rx^(rAlpha).*ry^(1 - rAlpha);

%Consider population of 50000 workers
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
%Parameterize the beta distribution
distCenter1     = [1,1];
distCenter2     = [1,1];


%Draw massPerPeriod workers and firms for each period such that there is a
%shift in the ability
adoptionProbW1   = betapdf(linspace(0,1,nx),distCenter1(1),distCenter1(2))';
adoptionProbF1   = betapdf(linspace(0,1,ny),distCenter1(1),distCenter1(2))';
adoptionProbW2   = betapdf(linspace(0,1,nx),distCenter2(1),distCenter2(2))';
adoptionProbF2   = betapdf(linspace(0,1,ny),distCenter2(1),distCenter2(2))';
%Do plot(adoptionProbW1) to see what the distributions are

%Sample the workers and firms in both periods. These become the workers and
%firms who are actually floating around in each period.
wName1          = find(sampleCDF(adoptionProbW1,wAdoptionNum));
wDist1          = workerDist(wName1);
wName2          = find(sampleCDF(adoptionProbW2,wAdoptionNum));
wDist2          = workerDist(wName2);
fName1          = find(sampleCDF(adoptionProbF1,fAdoptionNum));
fDist1          = firmDist(fName1);
fName2          = find(sampleCDF(adoptionProbF2,fAdoptionNum));
fDist2          = firmDist(fName2);

%Remember the bins stuff? We have 50 bins, but workers are from a
%continuous distribution. The averaging here screws things up a bit and
%further, screws things up across two periods.

%I think I can redo the
%code so that I can calculate the production function with a non-uniform
%density grid. 
%This would allow us to just sample the workers/firms from a fixed 50 point grid.
%I think it will do it, but take a look first before I implement that.
indProd1 = zeros(numBins,numBins);
indProd2 = zeros(numBins,numBins);
for i1 = 1:numBins
  ind1 = (i1-1)*workersPerBin + 1 : i1*workersPerBin;
  for i2 = 1:numBins
    ind2 = (i2-1)*firmsPerBin + 1 : i2*firmsPerBin;
    ind2 = vec(repmat(ind2,workersPerBin/firmsPerBin,1))';
    indProd1(i1,i2)        = mean(arrayfun(RP,wDist1(ind1),fDist1(ind2)));
    indProd2(i1,i2)        = mean(arrayfun(RP,wDist2(ind1),fDist2(ind2)));
  end
end

save

% This replication was performed on MATLAB R2012a
% This needs to be run from the \CODE folder
HLM('customProd',indProd1)
HLM('customProd',indProd2)
%Running this will create a folder Output where there are two saved files.
%Each of them plus stuff here will contain all the information we need to
%do our thing.

clear
load
load Output\customProd_customProd_000001_170512163349.mat
wTrueRank1    = RD.I.iNRRankAgg(:,1);
wEstRank1     = RD.I.iNRRankAgg(:,2);
wTrueBin1     = SimO.iNameX(wTrueRank1);
wEstBin1      = RD.I.iBin;
fTrueRank1    = vec(1:300);
fEstRank1     = RD.J.NROmega(:,2);
fTrueBin1     = ceil(vec(1:300)/6);
fEstBin1      = RD.J.jBin;
EconomyW1     = table(wName1,wTrueRank1,wEstRank1,wTrueBin1,wEstBin1);
EconomyF1     = table(fName1,fTrueRank1,fEstRank1,fTrueBin1,fEstBin1);
writetable(EconomyW1,'EconomyW1.csv')
writetable(EconomyF1,'EconomyF1.csv')

clear
load
load Output\customProd_customProd_000001_170512163843.mat
wTrueRank2    = RD.I.iNRRankAgg(:,1);
wEstRank2     = RD.I.iNRRankAgg(:,2);
wTrueBin2     = SimO.iNameX(wTrueRank2);
wEstBin2      = RD.I.iBin;
fTrueRank2    = vec(1:300);
fEstRank2     = RD.J.NROmega(:,2);
fTrueBin2     = ceil(vec(1:300)/6);
fEstBin2      = RD.J.jBin;
EconomyW2     = table(wName2,wTrueRank2,wEstRank2,wTrueBin2,wEstBin2);
EconomyF2     = table(fName2,fTrueRank2,fEstRank2,fTrueBin2,fEstBin2);
writetable(EconomyW2,'EconomyW2.csv')
writetable(EconomyF2,'EconomyF2.csv')

clear
addpath Source

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
distCenter1     = [2,5];
distCenter2     = [5,2];


%Draw massPerPeriod workers and firms for each period such that there is a
%shift in the ability
adoptionProbW1   = betapdf(linspace(0,1,nx),distCenter1(1),distCenter1(2))';
adoptionProbF1   = betapdf(linspace(0,1,ny),distCenter1(1),distCenter1(2))';
adoptionProbW2   = betapdf(linspace(0,1,nx),distCenter2(1),distCenter2(2))';
adoptionProbF2   = betapdf(linspace(0,1,ny),distCenter2(1),distCenter2(2))';
%Do plot(adoptionProbW1) to see what the distributions are

%Sample the workers and firms in both periods. These become the workers and
%firms who are actually floating around in each period.
wDist1          = workerDist(sampleCDF(adoptionProbW1,wAdoptionNum));
wDist2          = workerDist(sampleCDF(adoptionProbW2,wAdoptionNum));
fDist1          = firmDist(sampleCDF(adoptionProbF1,fAdoptionNum));
fDist2          = firmDist(sampleCDF(adoptionProbF2,fAdoptionNum));

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

% This replication was performed on MATLAB R2012a
% This needs to be run from the \CODE folder
HLM('customProd',indProd1)
HLM('customProd',indProd2)
%Running this will create a folder Output where there are two saved files.
%Each of them plus stuff here will contain all the information we need to
%do our thing.




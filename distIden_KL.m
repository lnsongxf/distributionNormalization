%Consider workers (1 to 15000) and firm IDs between 1 and 15000 of which 10k in each
%period
clear

%Consider RP (Real Production) defined on RX and RY 
rAlpha = 0.3;
RP = @(rx,ry)  0.5 + rx^(rAlpha).*ry^(1 - rAlpha);

%Consider population of 50000 workers
%30000 are in each one to match benchmark model
%There is a small issue with mixing workers and firms. Need to move
%carefully in understanding how selection affects the binning. Might need a
%theorem or two on continuity for indentification
nx              = 50000;
firmSize        = 100;
ny              = nx/firmSize;
rxmin           = 0;
rxmax           = 1;
rymin           = 0;
rymax           = 1;
wAdoptionNum    = 30000;
fAdoptionNum    = wAdoptionNum/firmSize;
distCenter1     = [2,5];
distCenter2     = [5,2];
numBins         = 50;
workersPerBin   = wAdoptionNum/numBins;
firmsPerBin     = fAdoptionNum/numBins;

%Distributions of the true underlying ability which is coming from some
%arbitrary distribution bound on [rxmin,rxmax] and [rymin,rymax]
workerDist = sort(rand(nx,1)*(rxmax - rxmin) + rxmin);
firmDist   = sort(rand(ny,1)*(rymax - rymin) + rymin);

%Draw massPerPeriod workers and firms for each period such that there is a
%shift in the ability
adoptionProbW1   = betapdf(linspace(0,1,nx),distCenter1(1),distCenter1(2))';
adoptionProbF1   = betapdf(linspace(0,1,ny),distCenter1(1),distCenter1(2))';
adoptionProbW2   = betapdf(linspace(0,1,nx),distCenter2(1),distCenter2(2))';
adoptionProbF2   = betapdf(linspace(0,1,ny),distCenter2(1),distCenter2(2))';

wDist1          = workerDist(sampleCDF(adoptionProbW1,wAdoptionNum));
wDist2          = workerDist(sampleCDF(adoptionProbW2,wAdoptionNum));
fDist1          = firmDist(sampleCDF(adoptionProbF1,fAdoptionNum));
fDist2          = firmDist(sampleCDF(adoptionProbF2,fAdoptionNum));

%Just start with a case where production function is just the appropriate
%average
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
addpath Source
HLM('customProd',indProd1)
HLM('customProd',indProd2)
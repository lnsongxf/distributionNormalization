%Consider workers (1 to 15000) and firm IDs between 1 and 15000 of which 10k in each
%period


%Consider RP (Real Production) defined on RX and RY 
rAlpha = 0.3;
RP = @(rx,ry)  rx^(rAlpha).*ry^(1 - rAlpha);

nxy             = 10000;
rxmin           = 0;
rxmax           = 1;
rymin           = 0;
rymax           = 1;
adoptionProb    = 0.7;
distCenter1     = [2,5];
distCenter2     = [5,2];
numBins         = 50;

%Distributions of the true underlying ability which is coming from some
%arbitrary distribution bound on [rxmin,rxmax] and [rymin,rymax]
workerDist = sort(rand(nxy,1)*(rxmax - rxmin) + rxmin);
firmDist   = sort(rand(nxy,1)*(rymax - rymin) + rymin);

%Draw massPerPeriod workers and firms for each period such that there is a
%shift in the ability
adoptionProb1   = betapdf(linspace(0,1,nxy),distCenter1(1),distCenter1(2))';
adoptionProb2   = betapdf(linspace(0,1,nxy),distCenter2(1),distCenter2(2))';

W1              = sampleCDF(adoptionProb1,adoptionProb*nxy);
wDist1          = workerDist(W1);
W2              = sampleCDF(adoptionProb2,adoptionProb*nxy);
wDist2          = workerDist(W2);
F1              = sampleCDF(adoptionProb1,adoptionProb*nxy);
fDist1          = firmDist(F1);
F2              = sampleCDF(adoptionProb2,adoptionProb*nxy);
fDist2          = firmDist(F2);

%Now, construct the production function based on the ranks of these
%workers. 
%Production function and hence worker wages are constructed at the
%individual level.
%Firms can be aggregates of jobs later on.
indProd1 = zeros(numel(wDist1),numel(fDist1));
indProd2 = zeros(numel(wDist2),numel(fDist2));
for i1 = 1:numel(wDist1)
  for i2 = 1:numel(fDist1)
    indProd1(i1,i2)        = RP(wDist1(i1),fDist1(i2));
    indProd2(i1,i2)        = RP(wDist2(i1),fDist2(i2));
  end
end

%Now construct the production function as averages of cells of firms of
%size 100


% This replication was performed on MATLAB R2012a
% This needs to be run from the \CODE folder
clear
addpath Source
HLM('ojs')

mkdir('Output\benchmark');
!move Output\benchmark* Output\benchmark\.
compileData('Output\benchmark')

makeProdError
benchmarkFigures
ojsFigures
matchqualityFigures
highbetaFigures
shortsampleFigures
smallfirmsFigures
plotDataFromIAB

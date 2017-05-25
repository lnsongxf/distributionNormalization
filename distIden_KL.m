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
RP{3} = @(x,y) (0.1 + (x-0.1+1)*y).*double(x<=0.1) + (0.1+((x-0.1)^2+y^2)^(1/2)).*double(x>0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Period 1 distribution normalized to a uniform
distCenter1     = [1,1];
%Period 2 distributions can be different from Period 1
distCenter2{1}     = [1,1];
distCenter2{2}     = [2,4];
distCenter2{3}     = [4,2];
save setup.mat



%This model

clear 
close all
TrueParam1                    = [1.7,1.3,0.4];
TrueParam2                    = [1.3,1.7,0.6];
TrueDist                      = dataset(vec(1:100000),sort(rand(100000,1)));
[~,~,TrueDist.Var3]           = unique(TrueDist.Var2,'sorted');
TrueDist.Properties.VarNames  = {'Name','Skill','Rank0'};
P1Sampling                    = betapdf(TrueDist.Skill,TrueParam1(1),TrueParam1(2));
P1Sampling                    = P1Sampling./max(P1Sampling);
Dist1                         = TrueDist(rand(100000,1)<=P1Sampling,:);
P1Sampling                    = betapdf(TrueDist.Skill,TrueParam2(1),TrueParam2(2));
P1Sampling                    = P1Sampling./max(P1Sampling);
Dist2                         = TrueDist(rand(100000,1)<=P1Sampling,:);
ProdFn1                       = @(x) x.^TrueParam1(3);
ProdFn2                       = @(x) x.^TrueParam2(3);
Dist1.Wage                    = ProdFn1(Dist1.Skill);
Dist2.Wage                    = ProdFn2(Dist2.Skill);
[~,~,Dist1.Rank1]             = unique(Dist1.Skill,'sorted');
[~,~,Dist2.Rank2]             = unique(Dist2.Skill,'sorted');
m1_TRUE                       = getMoments(Dist1.Wage);
m2_TRUE                       = getMoments(Dist2.Wage);

%Confirm that the distributions in the first and second half of the data is
%the Beta distributions used above.
[PHAT1, ~]                   = betafit(Dist1.Skill,0.05);
[PHAT2, ~]                   = betafit(Dist2.Skill,0.05);
PHAT1
TrueParam1
PHAT2
TrueParam2

%TRUE Decompositions.
m2UsingProd1_TRUE                 = getMoments(ProdFn1(Dist2.Skill));
m2UsingDist1_TRUE                 = getMoments(ProdFn2(Dist1.Skill));


%Nonparametric Approach
%Know that wages rank workers, and from there, simply construct production
%function on percentiles of workers. Consider the changes after normalizing
%based on ranks
%Estimate the production function nonparametrically in both periods.
nonParProd1            = getNonParProd(Dist1);
nonParProd2            = getNonParProd(Dist2);

%Estimate the inversion
%Obtain the rankings of workers in each half by wages
[~,~,Dist1.WageRank1]  = unique(Dist1.Wage,'sorted');
Dist1.WageRank1        = Dist1.WageRank1./max(Dist1.WageRank1);
[~,~,Dist2.WageRank2]  = unique(Dist2.Wage,'sorted');
Dist2.WageRank2        = Dist2.WageRank2./max(Dist2.WageRank2);
%Use only the overlaps because identification comes from there.
OverLap                = intersect(Dist1.Name,Dist2.Name);

%Construct the minimization problem.
F                      = @(x,xdata) kumaraswamyiCDF(xdata,x(1),x(2));
Fsumsquares            = @(x) sum((F(x,Dist2.WageRank2(ismember(Dist2.Name,OverLap))) - Dist1.WageRank1(ismember(Dist1.Name,OverLap))).^2);
% This is the inversion when applied to Overlapping W in period 2, yields
% ranking that is comparable to period 1.
xunc                   = fminunc(Fsumsquares,[0.5,0.5]);

% Alter the domains of the production function in period 2 as well as the
% distribution.
DomainProd1            = linspace(0,1,numel(nonParProd1));
DomainProd2            = linspace(0,1,numel(nonParProd2));
DomainProd2            = kumaraswamyiCDF(DomainProd2,xunc(1),xunc(2)); 
Dist2.WageRank2Norm    = kumaraswamyiCDF(Dist2.WageRank2,xunc(1),xunc(2)); 

m2UsingProd1_NonP      = getMoments(interp1(DomainProd1,nonParProd1,Dist2.WageRank2Norm));
m2UsingDist1_NonP      = getMoments(interp1(DomainProd2,nonParProd2,Dist1.WageRank1));

%Parametric Approach, No Specification Bias
%Econometrician correctly guesses that Wage function is Skill^Alpha
%And guesses that distribution is Beta.
%Wants to guess Beta distribution parameters and production function Alpha
options                = optimset('fmincon');
options.TolCon         = 1e-8;
options.TolFun         = 1e-8;
paramEst1              = fmincon(@(X) sum((simulateWages(X) - m1_TRUE).^2),TrueParam1,[],[],[],[],[0,0,0],[2,2,2],[],options);
paramEst2              = fmincon(@(X) sum((simulateWages(X) - m2_TRUE).^2),TrueParam2,[],[],[],[],[0,0,0],[2,2,2],[],options);

m2UsingProd1_Par       = getMoments(betarnd(paramEst2(1),paramEst2(2),100000,1).^paramEst1(3));
m2UsingDist1_Par       = getMoments(betarnd(paramEst1(1),paramEst1(2),100000,1).^paramEst2(3));

%Parametric Approach, Specification Bias in Production
%Econometrician correctly guesses that Wage function is Skill^Alpha
%And guesses that distribution is Beta.
%Wants to guess Beta distribution parameters and production function Alpha
paramEst1              = fmincon(@(X) sum((simulateWagesMisspec(X) - m1_TRUE).^2),TrueParam1,[],[],[],[],[0,0,0],[2,2,2],[],options);
paramEst2              = fmincon(@(X) sum((simulateWagesMisspec(X) - m2_TRUE).^2),TrueParam2,[],[],[],[],[0,0,0],[2,2,2],[],options);
m2UsingProd1_ParMS     = getMoments(betarnd(paramEst2(1),paramEst2(2),100000,1).^paramEst1(3));
m2UsingDist1_ParMS     = getMoments(betarnd(paramEst1(1),paramEst1(2),100000,1).^paramEst2(3));


disp('Moments from raw data.')
m1_TRUE
m2_TRUE
disp('Decompositions with true parameters.')
m2UsingProd1_TRUE                 = getMoments(ProdFn1(Dist2.Skill))
m2UsingDist1_TRUE                 = getMoments(ProdFn2(Dist1.Skill))
disp('Decompositions if we used the correct nonparametric normalization.')
m2UsingProd1_NonP      = getMoments(interp1(DomainProd1,nonParProd1,Dist2.WageRank2Norm))
m2UsingDist1_NonP      = getMoments(interp1(DomainProd2,nonParProd2,Dist1.WageRank1))
disp('Decompositions with correct parametric form.')
m2UsingProd1_Par       = getMoments(betarnd(paramEst2(1),paramEst2(2),100000,1).^paramEst1(3))
m2UsingDist1_Par       = getMoments(betarnd(paramEst1(1),paramEst1(2),100000,1).^paramEst2(3))
disp('Decompositions with incorrect parametric form.')
m2UsingProd1_ParMS     = getMoments(betarnd(paramEst2(1),paramEst2(2),100000,1).^paramEst1(3))
m2UsingDist1_ParMS     = getMoments(betarnd(paramEst1(1),paramEst1(2),100000,1).^paramEst2(3))

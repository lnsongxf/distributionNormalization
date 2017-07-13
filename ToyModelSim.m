function ToyModelSim()
  clear
  close all
  TrueParam1.Beta1              = 1.7;
  TrueParam1.Beta2              = 1.3;
  TrueParam1.ProdParam          = 0.4;
  TrueParam2.Beta1              = 1.3;
  TrueParam2.Beta2              = 1.7;
  TrueParam2.ProdParam          = 0.6;
  numWorkers                    = 100000;
  numGrid                       = 1000;
  measurementError              = 0.1;
  simRounds                     = 10;
  
  %Construct the "true" data
  %Name is worker id
  %Skill is the worker's skill in a quantifiable sense (eg:, IQ, schooling)
  %Rank is just an ordinal (eg: 2nd, 50th percentile)
  TrueDist                      = dataset(vec(1:numWorkers),sort(rand(numWorkers,1)));
  [~,~,TrueDist.Var3]           = unique(TrueDist.Var2,'sorted');
  TrueDist.Properties.VarNames  = {'Name','Skill','Rank0'};
  
  %We now select some workers to Period 1 and some workers to Period 2 so that
  %the distribution is Beta with different parameters
  Dist1       = getDist(TrueDist,TrueParam1);
  Dist2       = getDist(TrueDist,TrueParam2);
  
  disp('Verify that the distribution is indeed Beta')
  disp('True Period 1:')
  disp([TrueParam1.Beta1,TrueParam1.Beta2])
  disp('Estimated Period 1:')
  disp(betafit(Dist1.Skill,0.01))
  disp('True Period 2:')
  disp([TrueParam2.Beta1,TrueParam2.Beta2]);
  disp('Estimated Period 2:')
  disp(betafit(Dist2.Skill,0.01))
  
  %Construct production functions defined on cardinal skill
  ProdFn1                       = @(x) x.^TrueParam1.ProdParam;
  ProdFn2                       = @(x) x.^TrueParam2.ProdParam;
  
  %Wages are simply production plus some normal random error in logs.
  Dist1.Wage                    = exp(log(ProdFn1(Dist1.Skill)) + measurementError*randn(size(Dist1.Skill)));
  Dist2.Wage                    = exp(log(ProdFn2(Dist2.Skill)) + measurementError*randn(size(Dist2.Skill)));
  
  %In this model, wages rank workers, so we construct the ranking on [0,1]
  [~,~,Dist1.Rank]              = unique(Dist1.Wage,'sorted');
  Dist1.Rank                    = Dist1.Rank./max(Dist1.Rank);
  [~,~,Dist2.Rank]              = unique(Dist2.Wage,'sorted');
  Dist2.Rank                    = Dist2.Rank./max(Dist2.Rank);
  
  %Construct the 4 moments which we want to decompose.
  m1_TRUE                       = getMoments(Dist1.Wage);
  m2_TRUE                       = getMoments(Dist2.Wage);
  
  %The true counterfactuals are simply evaluating Prod2 using Skill1, and
  %Prod1 using Skill2.
  m2UsingProd1_TRUE                 = getMoments(ProdFn1(Dist2.Skill));
  m2UsingDist1_TRUE                 = getMoments(ProdFn2(Dist1.Skill));
  
  %Construct the minimization problem which is basically a rescaling of the
  %domain in Period 2 so that it is comparable to the domain in Period 1.
  %In otherwise, we are making the ranking of a worker in Period 2 as close
  %as possible to his ranking in Period 1.
  %This domain normalization will be used in both parametric and
  %non-parametric approaches.
  OverLap                = intersect(Dist1.Name,Dist2.Name); %Identification relies on overlapping workers.
  F                      = @(x,xdata) kumaraswamyiCDF(xdata,x(1),x(2));
  Fsumsquares            = @(x) sum((F(x,Dist2.Rank(ismember(Dist2.Name,OverLap))) - Dist1.Rank(ismember(Dist1.Name,OverLap))).^2);
  xunc                   = fminunc(Fsumsquares,[0.5,0.5]);
  
  %Nonparametric Approach
  %Estimate the production function nonparametrically in both periods.
  %This is simply a simple average of workers at every permille.
  nonParProd1            = getNonParProd(Dist1,numGrid);
  nonParProd2            = getNonParProd(Dist2,numGrid);
  
  % Alter domains of production function in period 2 and distribution.
  DomainProd1            = linspace(0,1,numel(nonParProd1));
  DomainProd2            = linspace(0,1,numel(nonParProd2));
  DomainProd2            = kumaraswamyiCDF(DomainProd2,xunc(1),xunc(2));
  Dist2.RankNorm         = kumaraswamyiCDF(Dist2.Rank,xunc(1),xunc(2));
  
  m2UsingProd1_NonP      = getMoments(interp1(DomainProd1,nonParProd1,Dist2.RankNorm));
  m2UsingDist1_NonP      = getMoments(interp1(DomainProd2,nonParProd2,Dist1.Rank));
  %Parametric Approach
  %We start the estimation at the TRUE value
  simMoments             = @(X) simulateWages(X,simRounds,numWorkers);
  paramEst1              = fminunc(@(X) sum((simMoments(X) - m1_TRUE).^2),[TrueParam1.Beta1,TrueParam1.Beta2,TrueParam1.ProdParam]);
  paramEst2              = fminunc(@(X) sum((simMoments(X) - m2_TRUE).^2),[TrueParam2.Beta1,TrueParam2.Beta2,TrueParam2.ProdParam]);
  
  m2UsingProd1_Par       = getMoments(betarnd(paramEst2(1),paramEst2(2),numWorkers,1).^paramEst1(3));
  m2UsingDist1_Par       = getMoments(betarnd(paramEst1(1),paramEst1(2),numWorkers,1).^paramEst2(3));
  
  %Parametric Approach, with Ranking as additional constraints
  %Start from the results of the parametric estimation, and simply alter the
  %domain of Period 2 distribution and production function to match Period 1
  %Covert the parametric estimation into estimation on percentiles
  for i1 = 1:100
    Idx                    = DomainProd1>((i1-1)*0.01) & DomainProd1<=(i1*0.01);
    Domain                 = DomainProd1(Idx);
    Prods                  = betainv(Domain,paramEst1(1),paramEst1(2));
    ParProd1(Idx,1)        = mean(Prods.^paramEst1(3));
    Prods                  = betainv(Domain,paramEst2(1),paramEst2(2));
    ParProd2(Idx,1)        = mean(Prods.^paramEst2(3));
  end
  m2UsingProd1_ParR      = getMoments(interp1(DomainProd1,ParProd1,Dist2.RankNorm));
  m2UsingDist1_ParR      = getMoments(interp1(DomainProd2,ParProd2,Dist1.Rank));
  
  disp('Moments from raw data.')
  m1_TRUE
  m2_TRUE
  disp('Decompositions with true parameters.')
  disp(sprintf('True              & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',             m2UsingProd1_TRUE(1),m2UsingProd1_TRUE(2),m2UsingProd1_TRUE(3),m2UsingProd1_TRUE(4)));
  disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',    m2UsingProd1_NonP(1),m2UsingProd1_NonP(2),m2UsingProd1_NonP(3),m2UsingProd1_NonP(4)));
  disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',       m2UsingProd1_Par(1),m2UsingProd1_Par(2),m2UsingProd1_Par(3),m2UsingProd1_Par(4)));
  disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',m2UsingProd1_ParR(1),m2UsingProd1_ParR(2),m2UsingProd1_ParR(3),m2UsingProd1_ParR(4)));
  disp(sprintf('True              & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',             m2UsingDist1_TRUE(1),m2UsingDist1_TRUE(2),m2UsingDist1_TRUE(3),m2UsingDist1_TRUE(4)));
  disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',    m2UsingDist1_NonP(1),m2UsingDist1_NonP(2),m2UsingDist1_NonP(3),m2UsingDist1_NonP(4)));
  disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',       m2UsingDist1_Par(1),m2UsingDist1_Par(2),m2UsingDist1_Par(3),m2UsingDist1_Par(4)));
  disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',m2UsingDist1_ParR(1),m2UsingDist1_ParR(2),m2UsingDist1_ParR(3),m2UsingDist1_ParR(4)));
  
end

function Dist  = getDist(TrueDist,TrueParam1)
  %Sample from the TrueDist(Uniform) so that we get Beta distribution
  P1Sampling   = betapdf(TrueDist.Skill,TrueParam1.Beta1,TrueParam1.Beta2);
  P1Sampling   = P1Sampling./max(P1Sampling);
  Dist        = TrueDist(rand(size(TrueDist.Name))<=P1Sampling,:);
end


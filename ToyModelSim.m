function ToyModelSim()
  clear
  close all
  TrueParam1.Beta1              = 1.7;
  TrueParam1.Beta2              = 1.3;
  TrueParam1.ProdParam          = 0.4;
  TrueParam2.Beta1              = 1.3;
  TrueParam2.Beta2              = 1.7;
  TrueParam2.ProdParam          = 0.6;
  numWorkers                    = 50000;
  numGrid                       = 1000;
  measurementError              = 0.1;
  simRounds                     = 10;
  randomStarts                  = 20;
  
  %Construct the "true" data
  %Name is worker id
  %Skill is the worker's skill in a quantifiable sense (eg:, IQ, schooling)
  %Rank is just an ordinal (eg: 2nd, 50th percentile)
  [TrueDist,Dist1,Dist2]        = constructData(numWorkers,TrueParam1,TrueParam2);
  
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
  Dist1.Wage                    = ProdFn1(Dist1.Skill) + measurementError*randn(size(Dist1.Skill));
  Dist2.Wage                    = ProdFn2(Dist2.Skill) + measurementError*randn(size(Dist2.Skill));
  
  %In this model, wages rank workers, so we construct the ranking on [0,1]
  [~,~,Dist1.Rank]              = unique(Dist1.Skill,'sorted');
  Dist1.Rank                    = Dist1.Rank./max(Dist1.Rank);
  [~,~,Dist2.Rank]              = unique(Dist2.Skill,'sorted');
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
  
  for i1 = 1:randomStarts
    startPoint = rand(1,6);
    startPoint = 1.2*[TrueParam1.Beta1 TrueParam1.Beta2 TrueParam1.ProdParam TrueParam2.Beta1 TrueParam2.Beta2 TrueParam2.ProdParam] + (startPoint-0.5);
    startPoint_vec(i1,:) = startPoint;
  end
  
  %Parametric Approach
  %COnsider outcomes of starting at randomly chosen points around the truth
  simMoments             = @(X) simulateWages(X,simRounds,numWorkers);
  
  options                           = optimoptions('fminunc','Algorithm','quasi-newton');
  options.Display                   = 'iter-detailed';
  options.StepTolerance             = 1e-6;
  options.FiniteDifferenceStepSize  = 1e-2;
  options.OptimalityTolerance       = 1e-6;
  
  for i1 = 1:randomStarts
    startPoint = startPoint_vec(i1,:);
    parFunc = @(X) sum((simMoments(X(1:3)) - m1_TRUE).^2) + ...
      sum((simMoments(X(4:6)) - m2_TRUE).^2);
    paramTogether = fminunc(parFunc, startPoint,options);
    paramEst1                      = paramTogether(1:3)
    paramEst2                      = paramTogether(4:6)
    paramEst1_vec(i1,:)            = paramTogether(1:3);
    paramEst2_vec(i1,:)            = paramTogether(4:6);
    m2UsingProd1_Par_vec(i1,:)     = getMoments(betarnd(paramEst2(1),paramEst2(2),numWorkers,1).^paramEst1(3));
    m2UsingDist1_Par_vec(i1,:)     = getMoments(betarnd(paramEst1(1),paramEst1(2),numWorkers,1).^paramEst2(3));
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Parametric Approach with overlapping worker moments
  % True Moment: (Dist1.Skill(I1)-Dist2.Skill(I2)) = 0
  % Sample Moment: betainv(Dist1.Rank(I1),X(1),X(2))-betainv(Dist2.Rank(I2),X(4),X(5)) = 0
  for i1 = 1:randomStarts
    startPoint = startPoint_vec(i1,:);
    parFunc = @(X) sum((simulateWages(X(1:3),simRounds,numWorkers) - m1_TRUE).^2) + ...
      sum((simulateWages(X(4:6),simRounds,numWorkers) - m2_TRUE).^2) + ...
      sum((getMoments(betainv(Dist1.Rank(ismember(Dist1.Name,OverLap)),X(1),X(2))) - ...
      getMoments(betainv(Dist2.Rank(ismember(Dist2.Name,OverLap)),X(4),X(5)))).^2);
    paramTogether = fminunc(parFunc,...
      startPoint,options);
    paramEst1_R                     = paramTogether(1:3)
    paramEst2_R                     = paramTogether(4:6)
    paramEst1R_vec(i1,:)            = paramTogether(1:3);
    paramEst2R_vec(i1,:)            = paramTogether(4:6);
    m2UsingProd1_ParR_vec(i1,:)     = getMoments(betarnd(paramEst2_R(1),paramEst2_R(2),numWorkers,1).^paramEst1_R(3));
    m2UsingDist1_ParR_vec(i1,:)     = getMoments(betarnd(paramEst1_R(1),paramEst1_R(2),numWorkers,1).^paramEst2_R(3));
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  m2UsingProd1_Par  = mean(m2UsingProd1_Par_vec);
  m2UsingDist1_Par  = mean(m2UsingDist1_Par_vec);
  m2UsingProd1_ParR = mean(m2UsingProd1_ParR_vec);
  m2UsingDist1_ParR = mean(m2UsingDist1_ParR_vec);
  
  save withgetMomentsBigSampleSpecBias.mat
  
%     disp('Decompositions with true parameters.')
%     disp(sprintf('True              & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',             m2UsingProd1_TRUE(1),m2UsingProd1_TRUE(2),m2UsingProd1_TRUE(3),m2UsingProd1_TRUE(4)));
%     disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',    m2UsingProd1_NonP(1),m2UsingProd1_NonP(2),m2UsingProd1_NonP(3),m2UsingProd1_NonP(4)));
%     disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',       m2UsingProd1_Par(1),m2UsingProd1_Par(2),m2UsingProd1_Par(3),m2UsingProd1_Par(4)));
%     disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',m2UsingProd1_ParR(1),m2UsingProd1_ParR(2),m2UsingProd1_ParR(3),m2UsingProd1_ParR(4)));
%     disp(sprintf('True              & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',             m2UsingDist1_TRUE(1),m2UsingDist1_TRUE(2),m2UsingDist1_TRUE(3),m2UsingDist1_TRUE(4)));
%     disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',    m2UsingDist1_NonP(1),m2UsingDist1_NonP(2),m2UsingDist1_NonP(3),m2UsingDist1_NonP(4)));
%     disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',       m2UsingDist1_Par(1),m2UsingDist1_Par(2),m2UsingDist1_Par(3),m2UsingDist1_Par(4)));
%     disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f & %7.4f \\\\ ',m2UsingDist1_ParR(1),m2UsingDist1_ParR(2),m2UsingDist1_ParR(3),m2UsingDist1_ParR(4)));
%   

%   disp('Decompositions with true parameters.')
%   disp(sprintf('True              & %7.4f & %7.4f  \\\\ ',             m2UsingProd1_TRUE(1),m2UsingProd1_TRUE(2)));
%   disp(sprintf('NonParametric     & %7.4f & %7.4f  \\\\ ',    m2UsingProd1_NonP(1),m2UsingProd1_NonP(2)));
%   disp(sprintf('Parametric        & %7.4f & %7.4f  \\\\ ',       m2UsingProd1_Par(1),m2UsingProd1_Par(2)));
%   disp(sprintf('Parametric + Rank & %7.4f & %7.4f  \\\\ ',m2UsingProd1_ParR(1),m2UsingProd1_ParR(2)));
%   disp(sprintf('True              & %7.4f & %7.4f  \\\\ ',             m2UsingDist1_TRUE(1),m2UsingDist1_TRUE(2)));
%   disp(sprintf('NonParametric     & %7.4f & %7.4f  \\\\ ',    m2UsingDist1_NonP(1),m2UsingDist1_NonP(2)));
%   disp(sprintf('Parametric        & %7.4f & %7.4f  \\\\ ',       m2UsingDist1_Par(1),m2UsingDist1_Par(2)));
%   disp(sprintf('Parametric + Rank & %7.4f & %7.4f  \\\\ ',m2UsingDist1_ParR(1),m2UsingDist1_ParR(2)));
%   
   disp('Decompositions with true parameters.')
    disp(sprintf('True              & %7.4f & %7.4f & %7.4f  \\\\ ',             m2UsingProd1_TRUE(1),m2UsingProd1_TRUE(2),m2UsingProd1_TRUE(3)));
    disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f  \\\\ ',    m2UsingProd1_NonP(1),m2UsingProd1_NonP(2),m2UsingProd1_NonP(3)));
    disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f  \\\\ ',       m2UsingProd1_Par(1),m2UsingProd1_Par(2),m2UsingProd1_Par(3)));
    disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f  \\\\ ',m2UsingProd1_ParR(1),m2UsingProd1_ParR(2),m2UsingProd1_ParR(3)));
    disp(sprintf('True              & %7.4f & %7.4f & %7.4f  \\\\ ',             m2UsingDist1_TRUE(1),m2UsingDist1_TRUE(2),m2UsingDist1_TRUE(3)));
    disp(sprintf('NonParametric     & %7.4f & %7.4f & %7.4f  \\\\ ',    m2UsingDist1_NonP(1),m2UsingDist1_NonP(2),m2UsingDist1_NonP(3)));
    disp(sprintf('Parametric        & %7.4f & %7.4f & %7.4f  \\\\ ',       m2UsingDist1_Par(1),m2UsingDist1_Par(2),m2UsingDist1_Par(3)));
    disp(sprintf('Parametric + Rank & %7.4f & %7.4f & %7.4f  \\\\ ',m2UsingDist1_ParR(1),m2UsingDist1_ParR(2),m2UsingDist1_ParR(3)));
  
end

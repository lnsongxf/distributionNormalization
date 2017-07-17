function [TrueDist,Dist1,Dist2] = constructData(numWorkers,TrueParam1,TrueParam2)
  TrueDist = dataset(vec(1:numWorkers),sort(rand(numWorkers,1)));
  [~,~,TrueDist.Var3]           = unique(TrueDist.Var2,'sorted');
  TrueDist.Properties.VarNames  = {'Name','Skill','Rank0'};
  
  %We now select some workers to Period 1 and some workers to Period 2 so that
  %the distribution is Beta with different parameters
  Dist1       = getDist(TrueDist,TrueParam1);
  Dist2       = getDist(TrueDist,TrueParam2);
end

function Dist  = getDist(TrueDist,TrueParam1)
  %Sample from the TrueDist(Uniform) so that we get Beta distribution
  P1Sampling   = betapdf(TrueDist.Skill,TrueParam1.Beta1,TrueParam1.Beta2);
  P1Sampling   = P1Sampling./max(P1Sampling);
  randNums     = rand(size(TrueDist.Name));
  randNums(1)  = 0;
  randNums(end)= 0;
  Dist         = TrueDist(randNums<=P1Sampling,:);
end
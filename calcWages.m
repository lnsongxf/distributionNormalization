function m = calcWages(X,Dist)
  Wages = betainv(Dist.Rank,X(1),X(2)).^X(3);
  m     = getMoments(Wages);
end
function m = simulateWages(X)
  rng(100)
  simWages  = betarnd(X(1),X(2),100000,1).^X(3);
  m         = getMoments(simWages);
end
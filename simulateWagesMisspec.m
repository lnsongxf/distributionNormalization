function m = simulateWagesMisspec(X)
  %Unused for now
  rng(100)
  simWages  = X(3)*betarnd(X(1),X(2),100000,1);
  m         = getMoments(simWages);
end
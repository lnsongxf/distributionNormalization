function m = simulateWages(X,simRound,NumAgents)
  m = zeros(simRound,4);
  for i1 = 1:simRound
    simWages  = betarnd(X(1),X(2),NumAgents,1).^X(3);
    m(i1,:)   = getMoments(simWages);
  end
  m = mean(m,1);
end
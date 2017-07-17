function m = getMoments(X)
  m(1) = mean(X);
  m(2) = moment(X,2);
  m(3) = moment(X,3);
end
function m = getMoments(X)
  m(1) = mean(X);
  m(2) = moment(X,2);
  m(3) = moment(X,3);
  m(4) = moment(X,4);
end
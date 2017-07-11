function Prod = getNonParProd(Dist)
  DistSort = sortrows(Dist,4);
  DistSort.RankWage = vec(1:size(DistSort,1));
  DistSort.RankWage = DistSort.RankWage./max(DistSort.RankWage)*100;
  DistSort.RankWage = ceil(DistSort.RankWage);
  Prod              = grpstats(DistSort.Wage,DistSort.RankWage,{'mean'});
end
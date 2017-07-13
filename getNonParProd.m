function Prod = getNonParProd(Dist,numGrid)
  DistSort          = sortrows(Dist,{'Rank'});
  DistSort.RankWage = vec(1:size(DistSort,1));
  DistSort.RankWage = DistSort.RankWage./max(DistSort.RankWage)*numGrid;
  DistSort.RankWage = ceil(DistSort.RankWage);
  Prod              = grpstats(DistSort.Wage,DistSort.RankWage,{'mean'});
end
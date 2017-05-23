%Recovering the original production function in Period 2
for iDist2 = 1:3
  for iProd2 = 1:3
    for iSet2 = 1:4
      addn     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
      bla = dir(['Output\*',addn,'.mat']);
      load(['Output\',bla(1).name]);
      load(['Output\',bla(2).name]);
      paramW = mapping{iDist2,iProd,iProd2,iSet,iSet2}.work_esti_bins.params;
      paramF = mapping{iDist2,iProd,iProd2,iSet,iSet2}.firm_esti_bins.params;
      %The objective is to uncover the underlying true production function
      %as much as possible. Underlying production plots the input
      %production function at the average productivity of workers and firms
      %in the bins
      TrueWProd = grpstats(wDist2,SimO.iNameX);
      TrueFProd = grpstats(fDist2,SimO.jNameY);
      %Now take RD.Y.Prod, and remap it to the orignal productivity
      EstWProd = kumaraswamyiCDF(linspace(0,1,50),paramW(1),paramW(2));
      EstFProd = kumaraswamyiCDF(linspace(0,1,50),paramF(1),paramF(2));
      figure
      surf(TrueFProd,TrueWProd,M.Prod);
      hold on
      surf(EstFProd,EstWProd,RD.Y.Prod,'FaceAlpha',0.5);
      surf(C.Grid,C.Grid,RD.Y.Prod,'FaceAlpha',0.2);
      hold off
      title([num2str(iDist2),num2str(iProd2),num2str(iSet2)])
    end
  end
end
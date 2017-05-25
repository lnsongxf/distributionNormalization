%Means, variances, correlations, parametric fit coefficients, differences
tPS.mean     = zeros(numel(mapping),1);
tPS.variance = zeros(numel(mapping),1);
tPS.corr     = zeros(numel(mapping),1);
tPS.diff     = zeros(numel(mapping),1);
tPS.maxdiff  = zeros(numel(mapping),1);
nPS          = tPS;
uPS          = tPS;
idxData      = zeros(numel(mapping),5);

icount = 0;
%Recovering the original production function in Period 2
for iProd1 = 1:3
  for iSet1 = 1:4
    for iDist2 = 1:3
      for iProd2 = 1:3
        for iSet2 = 1:4
          addn     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
          if exist(['.',filesep,'Output',filesep,addn,'.mat'],'file') > 0
            icount = icount + 1;
            idxData(icount,:) = [iProd1,iSet1,iDist2,iProd2,iSet2];
            bla = dir(['Output\*',addn,'.mat']);
            load(['Output\',bla(1).name]);
            load(['Output\',bla(2).name]);
            paramW = mapping{iDist2,iProd1,iProd2,iSet1,iSet2}.work_esti_bins.params;
            paramF = mapping{iDist2,iProd1,iProd2,iSet1,iSet2}.firm_esti_bins.params;
            %The objective is to uncover the underlying true production function
            %as much as possible. Underlying production plots the input
            %production function at the average productivity of workers and firms
            %in the bins
            TrueWProd = grpstats(wDist2,SimO.iNameX);
            TrueFProd = grpstats(fDist2,SimO.jNameY);
            %Now take RD.Y.Prod, and remap it to the orignal productivity
            EstWProd = kumaraswamyiCDF(linspace(0,1,numBins),paramW(1),paramW(2));
            EstFProd = kumaraswamyiCDF(linspace(0,1,numBins),paramF(1),paramF(2));
            EstWProd(1) = TrueWProd(1);
            EstFProd(1) = TrueFProd(1);
            EstWProd(end) = TrueWProd(end);
            EstFProd(end) = TrueFProd(end);
            
            %Compare this over the true distribution only.
            
            
            [mF,mW] = meshgrid(TrueWProd,TrueFProd);
            
            uProd = griddata(C.Grid,C.Grid,RD.Y.Prod,mF,mW);
            nProd = griddata(EstFProd,EstWProd,RD.Y.Prod,mF,mW);
            uProd(isnan(M.Prod)) = nan;
            nProd(isnan(M.Prod)) = nan;
            
            idx               = ~isnan(M.Prod) & ~isnan(nProd) & ~isnan(uProd);
            
            tPS.mean(icount,1) = nanmean(M.Prod(idx));
            nPS.mean(icount,1) = nanmean(nProd(idx));
            uPS.mean(icount,1) = nanmean(uProd(idx));
            
            tPS.variance(icount,1) = nanvar(M.Prod(idx));
            nPS.variance(icount,1) = nanvar(nProd(idx));
            uPS.variance(icount,1) = nanvar(uProd(idx));
            
            tPS.corr(icount,1)  = corr(M.Prod(idx),M.Prod(idx));
            nPS.corr(icount,1)  = corr(nProd(idx),M.Prod(idx));
            uPS.corr(icount,1)  = corr(uProd(idx),M.Prod(idx));
            
            tPS.diff(icount,1) = nanmean((M.Prod(idx) - M.Prod(idx)).^2);
            nPS.diff(icount,1) = nanmean((nProd(idx) - M.Prod(idx)).^2);
            uPS.diff(icount,1) = nanmean((uProd(idx) - M.Prod(idx)).^2);
            
            tPS.maxdiff(icount,1) = nanmax(abs(M.Prod(idx) - M.Prod(idx)))./tPS.mean(icount);
            nPS.maxdiff(icount,1) = nanmax(abs(nProd(idx) - M.Prod(idx)))./tPS.mean(icount);
            uPS.maxdiff(icount,1) = nanmax(abs(uProd(idx) - M.Prod(idx)))./tPS.mean(icount);
            
            
            %             if rand < 0.1
            %               figure
            %               plot(fT, [mF(idx),mW(idx)], M.Prod(idx))
            %               hold on
            %               plot(fU, [mF(idx),mW(idx)], nProd(idx))
            %               plot(fN, [mF(idx),mW(idx)], uProd(idx))
            %
            %
            %
            %               figure
            %               surf(TrueFProd,TrueWProd,M.Prod);
            %               hold on
            %               surf(TrueFProd,TrueWProd,nProd,'FaceAlpha',0.5);
            %               surf(TrueFProd,TrueWProd,uProd,'FaceAlpha',0.2);
            %               hold off
            %               title([num2str(iDist2),num2str(iProd2),num2str(iSet2)])
            %
            %               asdasd
            %             end
            
            
          end
        end
      end
    end
  end
end

save results.mat tPS nPS uPS idxData

% load results.mat 
% %Figures
% F1 = figure;
% set(F1,'Name','MeanDiff')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% plot(tPS.mean,nPS.mean,'kx')
% hold on
% plot(tPS.mean,uPS.mean,'ks')
% line([0.45 1.66],[0.45 1.66])
% axis([0.45 1.66 0.45 1.66])
% legend({'Normalized';'Uniform'},'Location','NorthWest')
% xlabel('Mean of True Production Function')
% ylabel('Mean of Estimated Production Function')
% title('Mean Differences of Estimated and True Production Function')
% 
% 
% F1 = figure;
% set(F1,'Name','VarDiff')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% plot(tPS.variance,nPS.variance,'kx')
% hold on
% plot(tPS.variance,uPS.variance,'ks')
% line([-0.1 0.1],[-0.1 0.1])
% axis([0 0.07 0 0.07])
% legend({'Normalized';'Uniform'},'Location','NorthWest')
% xlabel('Variance of True Production Function')
% ylabel('Variance of Estimated Production Function')
% title('Variance Differences of Estimated and True Production Function')
% 
% F1 = figure;
% set(F1,'Name','Corr1')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(nPS.corr)
% title('Correlations - Normalized')
% 
% F1 = figure;
% set(F1,'Name','Corr2')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(uPS.corr)
% title('Correlations - Uniform')
% 
% F1 = figure;
% set(F1,'Name','Norm1')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(nPS.diff)
% title('Euclidian Norm - Normalized')
% 
% F1 = figure;
% set(F1,'Name','Norm2')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(uPS.diff)
% title('Euclidian Norm - Uniform')
% 
% F1 = figure;
% set(F1,'Name','MaxDiff1')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(nPS.maxdiff)
% title('Max Diff - Normalized')
% 
% F1 = figure;
% set(F1,'Name','MaxDiff2')
% set(F1,'Units','Inches','outerposition',[0 0 11 11]);
% set(F1,'PaperOrientation','Landscape');
% set(F1,'PaperPositionMode','Auto');
% set(F1,'Clipping','off')
% hist(uPS.maxdiff)
% title('Max Diff - Uniform')



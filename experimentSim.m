iRun  = 1;
load setup.mat
%Consider all combinations of parameters to and parameters from.
% constructMapping
for iProd = 1:3
  for iSet = 1:8
    addn1     = ['Period1','_Prod',num2str(iProd),'_Dist1_Set',num2str(iSet)];
    if exist(['.',filesep,'Output',filesep,addn1,'.mat'],'file') > 0
      bla = dir(['Output\*',addn1,'.mat']);
      load(['Output\',bla(2).name]);
      P1Prod      = RD.Y.Prod;
      
      for iDist2 = 3:-1:1
        for iProd2 = 1:3
          for iSet2 = 1:8
            addn2     = ['Period2','_Prod',num2str(iProd2),'_Dist',num2str(iDist2),'_Set',num2str(iSet2)];
            if exist(['.',filesep,'Output',filesep,addn2,'.mat'],'file') > 0
              
              bla = dir(['Output\*',addn2,'.mat']);
              load(['Output\',bla(2).name]);
              P2Prod  = RD.Y.Prod;
              
              
              %What if only distribution had changed.
              %Take the production function estimated from the mapped distrib,
              %recompute using that.
              paramW    = mapping{iDist2,iProd,iProd2,iSet,iSet2}.work_esti_bins.params;
              paramF    = mapping{iDist2,iProd,iProd2,iSet,iSet2}.firm_esti_bins.params;
              
              %First period distribution is normalized
              TrueWProd     = grpstats(wDist1,SimO.iNameX);
              TrueFProd     = grpstats(fDist1,SimO.jNameY);
              [qF,qW]       = meshgrid(TrueFProd,TrueWProd);
              
              %Impost that the start and end is fine for stability
              EstWProd  = kumaraswamyiCDF(linspace(0,1,numBins),paramW(1),paramW(2));
              EstFProd  = kumaraswamyiCDF(linspace(0,1,numBins),paramF(1),paramF(2));
              EstWProd(1) = TrueWProd(1);
              EstFProd(1) = TrueFProd(1);
              EstWProd(end) = TrueWProd(end);
              EstFProd(end) = TrueFProd(end);
              if all(diff(EstWProd) > 0) == false
                error()
              end
              
              if all(diff(EstFProd) > 0) == false
                error()
              end
              
              %We need this interpolated on the productivities of period one.
              [meshF,meshW] = meshgrid(EstFProd,EstWProd);
              
              
              if iRun == 1
                %Use the counterfactual production function.
                idxgood   =~ isnan(P2Prod);
                indProd1  = griddata( meshF(idxgood),meshW(idxgood),P2Prod(idxgood),qF,qW);
                indProd1  = ndnanfilter(indProd1,@rectwin,[2,2]);
                indProd1(isnan(indProd1)) = 0;
                %Compute what happens if you use the second period production
                bla = ['P2Prod_t1','P',num2str(iProd),'S',num2str(iSet),'t2P',num2str(iProd2),'D',num2str(iDist2),'S',num2str(iSet)];
                try
                  getNLSInputs(indProd1,wAdoptionNum,bla,fAdoptionNum,numBins,iSet);
                end
              else
                %What if use second period distribution?
                bla = ['P2Dist_t1','P',num2str(iProd),'S',num2str(iSet),'t2P',num2str(iProd2),'D',num2str(iDist2),'S',num2str(iSet)];
                idxgood   =~ isnan(P1Prod);
                indProd2  = griddata(qF(idxgood),qW(idxgood),P1Prod(idxgood),meshF,meshW);
                indProd2  = ndnanfilter(indProd2,@rectwin,[2,2]);
                indProd2(isnan(indProd2)) = 0;
                try
                  getNLSInputs(indProd2,wAdoptionNum,bla,fAdoptionNum,numBins,iSet);
                end
              end
            end
          end
        end
      end
    end
  end
end
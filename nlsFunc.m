function results = nlsFunc(fEstBin1,fName2,wEstBin1,wName2,...
    fEstBin2,fTrueBin1,wEstBin2,wTrueBin1,...
    fEstRank1,fTrueBin2,wEstRank1,wTrueBin2,...
    fEstRank2,fTrueRank1,wEstRank2,wTrueRank1,...
    fName1,fTrueRank2,wName1,wTrueRank2,...
    numworkers,numfirms,numBins)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set Parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  aalpha       = 0.01;   % significance level for Kolmogorov-Smirnov Test
  x0           = [1 1];  % starting for optimiser
  
  for wf = {'work','firm'}
    for te = {'true','esti'}
      for rb = {'rank','bins'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xseq         = linspace(0,1,200);
        
        % Take overlapping workers/firms
        [~,IA,IB] = intersect(wName1,wName2);
        [~,JA,JB] = intersect(fName1,fName2);
        
        % Convert to percentiles
        if wf==1
          if rb==1
            if te==1
              datay = (wTrueRank1(IA))/numworkers;
              datax = (wTrueRank2(IB))/numworkers;
            else
              datay = (wEstRank1(IA))/numworkers;
              datax = (wEstRank2(IB))/numworkers;
            end
          else
            if te==1
              datay = (wTrueBin1(IA))/numBins;
              datax = (wTrueBin2(IB))/numBins;
            else
              datay = (wEstBin1(IA))/numBins;
              datax = (wEstBin2(IB))/numBins;
            end
          end
        else
          if rb==1
            if te==1
              datay = (fTrueRank1(JA))/numfirms;
              datax = (fTrueRank2(JB))/numfirms;
            else
              datay = (fEstRank1(JA))/numfirms;
              datax = (fEstRank2(JB))/numfirms;
            end
          else
            if te==1
              datay = (fTrueBin1(JA))/numBins;
              datax = (fTrueBin2(JB))/numBins;
            else
              datay = (fEstBin1(JA))/numBins;
              datax = (fEstBin2(JB))/numBins;
            end
          end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parametric estimation of CDF of x for Period 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Set up nonlinear least squares loss function
        F = @(x,xdata) kumaraswamyiCDF(xdata,x(1),x(2));
        Fsumsquares = @(x) sum((F(x,datax) - datay).^2);
        
        % Solve for parameters of inverse CDF of x for Period 2
        xunc = fminunc(Fsumsquares,x0);
        % display(xunc)
        
        est_cdf      = kumaraswamyCDF(xseq,xunc(1),xunc(2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        % figure;
        % hold on;
        % plot(xseq,true_cdf);
        % plot(xseq,est_cdf,'r');
        % title('CDF True v. Estimated');
        % legend('True','Estimated','Location','southeast');
        % hold off;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kolmogorov-Smirnov Test for equality in distributions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [h,pval] = kstest2(true_cdf,est_cdf,'alpha',aalpha);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spec = [wf,'_',te,'_',rb];
        results.(spec).h    = h;
        results.(spec).pval = pval;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Screen display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % disp('Kolmogorov-Smirnov Test for equality in distributions')
        % if h==1
        %     disp(['Reject null of equality at ',num2str(100*aalpha),'% level'])
        % else
        %     disp(['   Fail to reject null of equality at ',num2str(100*aalpha),'% level'])
        % end
        % disp(['   p-value: ',num2str(pval)]);
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
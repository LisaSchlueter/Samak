% old stuff in PlotFit for ring option
  % get right of horizontal line in errorbars
            if numel(hebar)==1
                hebar.CapSize = 0;
            else
                c = jet(numel(hebar));
                for i=1:numel(hebar)
                    hebar(i).CapSize = 0;
                    lfit(i).Color = c(i,:);
                end
            end
            
            
             else %rings
                myleg.FontSize = LocalFontSize;
                ringleg = split(sprintf('ring:%.0f ',obj.RingList),' ');
                switch FitResultsFlag
                    case 'ON'
                        if contains(obj.fixPar,' 1 ')
                            mnuleg = sprintf('m_\\beta fixed,    ');
                        else
                            mnuleg = sprintf('m^2_\\beta %.2f \t\\pm %.2g eV^2,    ',obj.FitResult.par(1)+obj.ModelObj.mnuSq_i,obj.FitResult.err(1));
                        end
                        chi2leg = sprintf('\\chi2 = %.1f/%.0f dof \n',obj.FitResult.chi2min,obj.FitResult.dof);
                        e0leg = sprintf('E_0 = %.2f \t\\pm %.2g eV\n',obj.FitResult.par(2)+obj.ModelObj.Q_i,obj.FitResult.err(2));
                        fitlegCommon = [chi2leg,mnuleg,e0leg];
                        dataleg = [sprintf('data \n'),fitlegCommon];
                        norms_fit = obj.FitResult.par(3+obj.ModelObj.nPixels:3+2*obj.ModelObj.nPixels-1);
                        norms_fit_err = obj.FitResult.err(3+obj.ModelObj.nPixels:3+2*obj.ModelObj.nPixels-1);
                        qUoffset_fit     = obj.FitResult.par(2*obj.ModelObj.nPixels+9:3*obj.ModelObj.nPixels+8);
                        qUoffset_fit_err = obj.FitResult.err(2*obj.ModelObj.nPixels+9:3*obj.ModelObj.nPixels+8);
                        % Background (pixelwise)
                        bcks_fit     = obj.FitResult.par(3:2+obj.ModelObj.nPixels);
                        bcks_fit_err = obj.FitResult.err(3:2+obj.ModelObj.nPixels);
                        bnleg = string(sprintf('(ring %.0f): B = %.1f \\pm %.1f mcps,     N = %.2f \\pm %.2f ,     \\Delta qU = %.2f \\pm %.2g eV',...
                            [obj.ModelObj.FPD_RingList; (obj.ModelObj.BKG_RateSec_i + bcks_fit)*1e3;bcks_fit_err*1e3;norms_fit+1;norms_fit_err;qUoffset_fit;qUoffset_fit_err]));
                        fitleg =split(bnleg,'eV');
                        myleg = legend([hdata(1), lfit'],dataleg,fitleg{1:obj.nRings},'Location','Northeast') ;
                        myleg.FontSize = LocalFontSize;
                      
                    case 'OFF'
                        myleg = legend([hdata(1), lfit'],'data',ringleg{1:obj.nRings},'Location','Northeast') ;
                        myleg.FontSize = LocalFontSize;
                end
            end
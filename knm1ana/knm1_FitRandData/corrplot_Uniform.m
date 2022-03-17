function corrplot_Uniform(FitPar)
% correlation plot for 4 fit parameter (Uniform fit)

f1 = figure('Units','normalized','Position',[-1 0.1 0.7 0.8]);

% plot everything with respect to mean
FitPar(1,:) = FitPar(1,:) - mean(FitPar(1,:));
FitPar(2,:) = FitPar(2,:) - mean(FitPar(2,:));
FitPar(3,:) = FitPar(3,:) - mean(FitPar(3,:));
FitPar(4,:) = (FitPar(4,:) - mean(FitPar(4,:))).*100;

CorrMat = corrcoef(FitPar');
[~,StdCorrMat,~] = corrcoeff_err(FitPar);
 %% plot
for i=1:4
    for j=1:4
        if i==j
            s1 = subplot(4,4,4*(i-1)+j);
            histogram(FitPar(i,:),'EdgeColor','none','FaceAlpha',1,'FaceColor',rgb('DeepSkyBlue'))
            MakePretty
            yticks([]);
             ax = gca;    
        elseif i>j
            s1 = subplot(4,4,4*(i-1)+j);
            dscatter(FitPar(j,:)',FitPar(i,:)')
            ax = gca;
            MakePretty
            hold on;
            pnone = plot(NaN,NaN,'wo','MarkerFaceColor','none','MarkerEdgeColor','none');
            if round(StdCorrMat(i,j),2)==0
                leg = legend(pnone,sprintf('%.2f \\pm %.3f',CorrMat(i,j),StdCorrMat(i,j)),'Location','northwest');
            else
                leg = legend(pnone,sprintf('%.2f \\pm %.2f',CorrMat(i,j),StdCorrMat(i,j)),'Location','northwest');
            end
            PrettyLegendFormat(leg);
            leg.ItemTokenSize = [0,0];
            leg.TextColor = rgb('DeepPink');
            
        else
            continue
        end
        
        
        %% x-y labels
        if i==4 % bottom row
            if j==1
                lStr = sprintf('\\Delta{\\itm}_\\nu^{ 2} (eV^{ 2})');
                ylabel(sprintf('\\Delta{\\itN}_{sig.} (%%)'));
                if max(FitPar(1,:))<2.2 %KNm2
                  ax.YLabel.Position(1) = -2.3;
                else
                     ax.YLabel.Position(1) = -5.4;
                end
            elseif j==2
                lStr = sprintf('\\Delta{\\itE}_0^{ fit} (eV)');
            elseif j==4
                lStr = sprintf('\\Delta{\\itN}_{sig.} (%%)');
            elseif j==3
                lStr = sprintf('\\Delta{\\itB}_{base} (mcps)');
            end
            xlabel(lStr); 
        elseif j==1 % first column
            if i==1
                lStr = sprintf('\\Delta{\\itm}_\\nu^{ 2} (eV^{ 2})');
            elseif i==2
                lStr = sprintf('\\Delta{\\itE}_0^{ fit} (eV)');
            elseif i==4
                lStr = sprintf('\\Delta{\\itN}_{sig.} (%%)');
            elseif i==3
                lStr = sprintf('\\Delta{\\itB}_{base} (mcps)');
            end
            ylabel(lStr);
            if max(FitPar(1,:))>2.2
                % KNM1
               % if i==1
                    ax.YLabel.Position(1) = -5.4;
              %  elseif i==3
               %     ax.YLabel.Position(1) = -4.92;
              %  end
            else
                % KNM2
                if i==1
                    ax.YLabel.Position(1) = -2.3;
                elseif i==3
                    ax.YLabel.Position(1) = -2.3;
                end
            end
        end
        
        %% ticks
        if j==2  % endpoint points (2nd column)
            ax.XAxis.Exponent = 0;
        elseif i==2 % endpoint points (2nd row)
            ax.YAxis.Exponent = 0;
        end
        
        if i<4 % everything, but bottom row
            xticks([]);
        end
        
        %% move panels closer  together
        if j>1 % everything, but first column
            yticks([]);
            if j==2
            ax.Position(1) = 0.288;
            elseif j==3
                 ax.Position(1) = 0.446;
            elseif j==4
                 ax.Position(1) = 0.603;
            end
        end
        
        % axis limits
        if max(FitPar(1,:))>2.2
            % KNM1
            if j==1
                xlim([-3.5 3.1])
            elseif j==2
                xlim([-0.25 0.25])
                xticks([-0.2 0 0.2]);
            elseif j==4
                xlim([-1.3 1.5])
            elseif j==3
                xlim([-3 3])
            end
            
            if i==1 && j~=i
                ylim([-3.5 3.1])
            elseif i==2 && j~=i
                ylim([-0.25 0.25])
                yticks([-0.2 0 0.2]);
            elseif i==4 && j~=i
                ylim([-1.3 1.5])
            elseif i==3&& j~=i
                ylim([-3 3])
            end
        else   
            % KNM2
            if j==1
                xlim([-1.5 1.5])
            elseif j==2
                xlim([-0.14 0.14])
                xticks([-0.1 0 0.1]);
            elseif j==4
                xlim([-1.1 1.1])
                xticks([-0.7 0 0.7])
            elseif j==3
                xlim([-2.5 2.5])
            end
            
            if i==1 && j~=i
                ylim([-1.5 1,5])
            elseif i==2 && j~=i
                ylim([-0.14 0.14])
                xticks([-0.1 0 0.1]);
            elseif i==4 && j~=i
                ylim([-1.1 1.1])   
                yticks([-0.7 0 0.7])
                if j~=1
                    yticklabels('');
                end
            elseif i==3&& j~=i
                ylim([-2.5 2.5])
            end
        
        end
    end
end
end


function MakePretty
PrettyFigureFormat('FontSize',22);
ax = gca;
ax.Position(3) = 0.15;
ax.Position(4) = 0.21;
end


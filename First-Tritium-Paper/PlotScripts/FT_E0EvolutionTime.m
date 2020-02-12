% Fit every run in RunList
Mode = 'Plot';
myRuns = 'FTpaper';
mychi2 = 'chi2Stat';
exclDataStart = 13;
if exclDataStart==1
    belowE0 = 1600;
elseif exclDataStart==7
    belowE0 = 400;
elseif exclDataStart==9
    belowE0 = 200;
elseif exclDataStart==13
    belowE0 = 100;
end
%%
save_name = [getenv('SamakPath'),sprintf('studies/local_LisaThesisPlots/results/FitRunList_%s_%.0feV_%s.mat',mychi2,belowE0,myRuns)];

switch Mode
    case 'Compute'
        MRA = MultiRunAnalysis('RunList',myRuns,'chi2',mychi2,'fixPar','1 5 6 7 8 9 10 11',...
            'DataEffCor','RunSummary',...
            'ELossFlag','Abdurashitov','exclDataStart',exclDataStart,'FSDFlag','SAENZ','DataType','Real',...
             'NonPoissonScaleFactor',1);
        try
            FitResults = MRA.FitRunList('RecomputeFlag','ON');
        catch
            FitResults = MRA.FitRunList('RecomputeFlag','ON');
        end
        save(save_name,'FitResults');
    case 'Plot'
        load(save_name);
        M = MultiRunAnalysis('RunList',myRuns,'chi2',mychi2,'exclDataStart',exclDataStart,...
            'ELossFlag','Abdurashitov','FSDFlag','SAENZ','DataType','Real','fixPar','1 5 6 7 8 9 10 11',...
            'NonPoissonScaleFactor',1);
        %Values per Run
        Par = 'E0';
        E0 = FitResults.(Par); %E0 = E0(M.SingleRunData.Select_all);
        E0Err = FitResults.([Par,'Err']); %E0Err = E0Err(M.SingleRunData.Select_all);
        if strcmp(Par,'B')
            E0 = E0*1e3;
            E0Err = E0Err*1e3;
        elseif strcmp(Par,'N')
            E0 = E0+1;
        end
        % Weighted mean and error of mean
        switch mychi2
            case 'chi2Stat'
                meanE0      = wmean(E0,1./E0Err.^2);
                ErrOfMeanE0 = sqrt(1/sum(1./E0Err.^2));
                chi2min = sum((E0-meanE0).^2./E0Err.^2);
                pvalue  = 1-chi2cdf(chi2min,numel(E0)-1);
                fprintf(2,'%s \nchi2 = %.1f (%.0f dof)\n p-value = %.2f \n',Par, chi2min, numel(E0)-1,pvalue);
            case {'chi2CM','chi2CMShape'}
                CorrCoeff = (0:0.1:1);
                meanE0_v = zeros(numel(CorrCoeff,1));
                ErrOfMeanE0_v = zeros(numel(CorrCoeff,1));
                mychi2min = zeros(numel(CorrCoeff,1));
                mydof = zeros(numel(CorrCoeff,1));
                for i=1:numel(CorrCoeff)
                myCorrCoeff = 1;%10;
                %i=myCorrCoeff;
                try
                    [meanE0_v(i), ErrOfMeanE0_v(i), mychi2min(i), mydof(i), SysCM, ResultCM] = ...
                        M.FitRunList_AveragePar('ReFit','OFF','CorrCoeff',CorrCoeff(i),'Parameter',Par);
                catch
                    [meanE0_v(i), ErrOfMeanE0_v(i), mychi2min(i), mydof(i), SysCM, ResultCM] = ...
                        M.FitRunList_AveragePar('ReFit','ON','CorrCoeff',CorrCoeff(i),'Parameter',Par);
                end
                end
               
                meanE0 = meanE0_v(myCorrCoeff);
                ErrOfMeanE0 = ErrOfMeanE0_v(myCorrCoeff);
                if strcmp(Par,'B')
                    meanE0 = meanE0*1e3;
                    meanE0_v = meanE0_v*1e3;
                    ErrOfMeanE0 = ErrOfMeanE0*1e3;
                    ErrOfMeanE0_v = ErrOfMeanE0_v*1e3;
                elseif strcmp(Par,'N')
                    meanE0 = meanE0+1;
                    meanE0_v = meanE0_v+1;
                end
                close
                 fprintf(2,'%s \nchi2 = %.1f (%.0f dof)\n p-value = %.2f \n',...
                     Par, mychi2min(myCorrCoeff), numel(E0)-1,1-chi2cdf(mychi2min(myCorrCoeff),numel(E0)-1));
        end
        
        nRuns = numel(E0);
        CumTime = cumsum(M.SingleRunData.TimeSec(M.SingleRunData.Select_all)./(60*60));%in hours
        CumTime = CumTime([1,round(nRuns/4),round(nRuns/2),round(nRuns*3/4),nRuns]);
        %%
        fig1 = figure('Renderer','painters');
        set(fig1,'units','centimeters','pos',[0.1, 0.1,17.3,5.5]);%17.3,5.5
        MultiMean = 'OFF';
        
        x = hours(M.SingleRunData.StartTimeStamp-M.SingleRunData.StartTimeStamp(1));
        xIndex = x<150;
        x1 = x(x<150);
        x2 = x(x>150);
        xdiff = min(x2)-max(x1)-10;
        x2 = x2-xdiff;
        
        if strcmp(mychi2,'chi2Stat') || strcmp(MultiMean,'OFF')
            l = plot(linspace(min(x)-2,max(x)+2,10),meanE0.*zeros(1,10),...
                '-','LineWidth',1,'Color',rgb('Black')); hold on
        elseif strcmp(mychi2,'chi2CM') || strcmp(mychi2,'chi2CMShape') || strcmp(MultiMean,'ON')
            l=  plot(0:nRuns+1,meanE0_v.*ones(nRuns+2,1),...
                '-','LineWidth',1); hold on
            set(l,{'color'},num2cell(jet(numel(l)),2));
            l = l(end);
        end
        
        e1 = errorbar(x1,E0(xIndex)-meanE0,E0Err(xIndex),...
            'o','MarkerSize',3,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1,'Color',rgb('DarkCyan'));
        hold on;
        e1.CapSize = 0;
        e1 = errorbar(x2,E0(~xIndex)-meanE0,E0Err(~xIndex),...
            'o','MarkerSize',3,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1,'Color',rgb('DarkCyan'));
        %p = plot(max(x1).*ones(2,1)+5,[-10 10],'-','Color',rgb('Silver'),'LineWidth',6);
        
        %         e1 = errorbar(1:nRuns,E0-meanE0,E0Err,...
        %             'o','MarkerSize',10,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',2,'Color',rgb('DarkCyan'));
        ylim([-8,7]);
        e1.CapSize = 0;
        myxlim=get(gca,'XLim');
        myylim=get(gca,'YLim');
        Pos_tmp = get(gcf,'Position');
        if Pos_tmp(3)>10
            ht = text(0.49*myxlim(1)+0.5*myxlim(2),0.96.*myylim(1),sprintf('\\mid\\mid'));
        else
            ht = text(0.4*myxlim(1)+0.405*myxlim(2),0.96.*myylim(1),sprintf('\\mid\\mid'));
        end
        set(ht,'Rotation',-25)
        set(ht,'FontSize',10)
        set(ht,'FontWeight','bold');
        set(ht,'FontName','Arial');
        
        if  strcmp(mychi2,'chi2Stat')
            meanleg = sprintf('< {\\itE}_0^{fit} > = %.2f eV , std = %.2f eV \n\\chi2 (dof) = %.1f (%.0f) \\rightarrow p-value = %.2f',meanE0,std(E0),chi2min, numel(E0)-1,pvalue);
             leg = legend([e1,l],'Scanwise fits (statistics only)','Weighted mean','Location','southwest');  
        elseif strcmp(MultiMean,'OFF')
             leg = legend([e1,l],'runwise fit',sprintf('weighted mean \\rho = %.1f ',CorrCoeff(myCorrCoeff)),'Location','southwest');    
        elseif  strcmp(mychi2,'chi2CM') || strcmp(mychi2,'chi2CMShape') || strcmp(MultiMean,'ON')
             leg = legend([e1,l],'runwise fit',sprintf('weighted mean \\rho = %.1f - %.1f',min(CorrCoeff),max(CorrCoeff)),'Location','southwest');
        end
        legend boxoff
        
        xlabelMode = 3;
        switch xlabelMode
            case 1
                xticks([1,round(nRuns/4),round(nRuns/2),round(nRuns*3/4),nRuns]);
                xticklabels(string(round(CumTime,0)));
                xlabel('Time (hours)')
            case 2
                M.SingleRunData.StartTimeStamp.Format = 'dd/MM/yy';
                t = M.SingleRunData.StartTimeStamp;
                xticks([1,round(nRuns/8,0),round(nRuns/4),round(nRuns*3/8),round(nRuns/2),round(nRuns*5/8),round(nRuns*3/4),round(nRuns*7/8),nRuns]);
                xticklabels(string(t(xticks)));
                xtickangle(25);
                xlabel('Time');
            case 3
             xlim([min(x1)-4,max(x2)+4])  
             xticksIndex = xticks<=max(x1);
             myxticks = xticks;
             str1 = string(myxticks(xticksIndex));
             str2 = string(round(myxticks(~xticksIndex)+xdiff,0));
             xticklabels([str1,str2]);
             xlabel('Time (hours)');
        end
        
        FTpaperFormat;
        set(gca,'FontSize',9);
        set(get(gca,'XLabel'),'FontSize',9);
        set(get(gca,'YLabel'),'FontSize',9);

        switch Par
            case 'E0'
                 ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle  (eV)'));% - < E_{0_{eff.}} > 
                 yticklabels(string(round(get(gca,'YTick'),5)));
                 leg.Location = 'southwest';
                % leg.NumColumns = ;
            case 'B'
                ylabel('B (mcps)');
                leg.Location = 'northeast';
            case 'N'
                ylabel('N');
        end
       
        leg.FontSize = get(gca,'FontSize');
      
        fig_str = sprintf('FitRunListPlot_%s_%s_%.0feV_%.0fruns_xlabel%.0f',Par,mychi2,belowE0,numel(M.StackedRuns),xlabelMode);
        export_fig(fig1,[getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/plots/',fig_str,'.pdf']);
       % print(fig1,[getenv('SamakPath'),'studies/local_LisaThesisPlots/plots/png/',fig_str,'.png'],'-dpng');
        
        %%
%         chi2min = FitResults.chi2min; dof = FitResults.dof(1);
%         chi2plot = linspace(min(chi2min)-5,max(chi2min)+5,100);
%         chi2dist = chi2pdf(chi2plot,dof);
%         
%         fig2 = figure('Renderer','opengl');
%          set(fig2,'units','normalized','pos',[0.1, 0.1,0.5,0.5]);
%          
%          h1 = histogram(chi2min/dof,6);
%          h1.FaceColor = rgb('CadetBlue'); h1.LineWidth=1;
%          hold on;
%          plot(chi2plot/dof,(chi2dist.*h1.BinWidth.*nRuns*dof),'LineWidth',4,'Color',rgb('GoldenRod'))
%          PrettyFigureFormat;
%          xlabel(sprintf('red. \\chi^2(%.0f dof)',dof));
%          ylabel('runs');
%          set(gca,'FontSize',20);
%         xlim([0 2.5]);
%          
%         fig_str = sprintf('FitRunListChi2_%s_%.0feV_70runs',mychi2,belowE0);
%         publish_figurePDF(fig2,['../local_LisaThesisPlots/plots/pdf/',fig_str,'.pdf']);
%         print(fig2,['../local_LisaThesisPlots/plots/png/',fig_str,'.png'],'-dpng');
end
% %% Plot E0 as a function of CorrCoeff
% M.chi2 = 'chi2CMShape';
% M.Fit;
% E0Stacked = M.FitResult.par(2)+M.ModelObj.Q_i;
% E0StackedErr = M.FitResult.err(2);
% save_stat = strrep(save_name,'chi2CMShape','chi2Stat');
% d = importdata(save_stat);
% M.chi2 = 'chi2Stat';
% M.Fit;
% E0StackedStat = M.FitResult.par(2)+M.ModelObj.Q_i;
% E0StackedErrStat = M.FitResult.err(2);
%%
Par = 'E0';
if strcmp(Par,'E0') && strcmp(mychi2,'chi2CMShape')
      fig12 = figure('Renderer','painters'); %Endpoint
        set(fig12, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
    [l,a] = boundedline(CorrCoeff,repmat(E0Stacked,numel(CorrCoeff),1)-E0Stacked,repmat(E0StackedErr,numel(CorrCoeff),1));
    l.LineWidth = 3; l.Color = rgb('GoldenRod');
    a.FaceColor = rgb('GoldenRod');
    a.FaceAlpha = 0.5;
    hold on;
  % [lstat, astat] = boundedline(CorrCoeff,repmat(E0StackedStat,1,numel(CorrCoeff))-E0Stacked,repmat(E0StackedErrStat,1,numel(CorrCoeff)));
  %  lstat.LineWidth = 3; lstat.Color = rgb('GoldenRod');
   % astat.FaceColor = rgb('GoldenRod');
   % astat.FaceAlpha = 0.5;
  %  hold on;
   e=  errorbar(CorrCoeff,meanE0_v-E0Stacked,ErrOfMeanE0_v,'--o',...
       'MarkerSize',10,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',2,'Color',rgb('DarkCyan'));
   % eStat=  errorbar(CorrCoeff,repmat(wmean(d.E0,d.E0Err),1,numel(CorrCoeff))-E0Stacked,repmat(sqrt(1/sum(1./d.E0Err.^2)),1,numel(CorrCoeff)),'--o',...
     %  'MarkerSize',10,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2,'Color',rgb('FireBrick'));
   
   PrettyFigureFormat;
   xlabel(sprintf('correlation coeffcient \\rho'))
   ylabel('\Delta E_{0_{eff.}} (eV)');
   leg = legend([l,e],'stacked','average','Location','northwest');
   legend boxoff;
   leg.FontSize = 20;
   set(gca,'FontSize',20);
   ylim([-0.6, 0.35]);
   yticklabels(string(round(get(gca,'YTick'),5)));
    
   fig_str = sprintf('ComparisonRunCombi_StackedAverage_%s_%s_%.0feV',mychi2,myRuns,belowE0);
   publish_figurePDF(fig12,[getenv('SamakPath'),'studies/local_LisaThesisPlots/plots/pdf/',fig_str,'.pdf']);
   print(fig12,[getenv('SamakPath'),'studies/local_LisaThesisPlots/plots/png/',fig_str,'.png'],'-dpng','-r400');
end
%% Plot E0 distribution


fig2 = figure('Renderer','painters');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.5,0.5]);

h1 = histfit(E0,20,'norm');
h1.FaceColor = rgb('CadetBlue'); h1.LineWidth=1;
hold on;
plot(chi2plot/dof,(chi2dist.*h1.BinWidth.*nRuns*dof),'LineWidth',4,'Color',rgb('GoldenRod'))
PrettyFigureFormat;
xlabel(sprintf('red. \\chi^2(%.0f dof)',dof));
ylabel('runs');
set(gca,'FontSize',20);
xlim([0 2.5]);

fig_str = sprintf('FitRunListChi2_%s_%.0feV_70runs',mychi2,belowE0);
publish_figurePDF(fig2,[getenv('SamakPath'),'studies/local_LisaThesisPlots/plots/pdf/',fig_str,'.pdf']);
print(fig2,[getenv('SamakPath'),'studies/local_LisaThesisPlots/plots/png/',fig_str,'.png'],'-dpng');
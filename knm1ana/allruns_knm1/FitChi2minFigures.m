%
% Fit Endpoint Stability with a Constant / Constant + Slope
% T. Lasserre
% Last Modified, June 3 2019
%

% Define Data Set
t          =  (datenum(Real.SingleRunData.StartTimeStamp)-datenum(Real.SingleRunData.StartTimeStamp(1)));%days
chi2min    =  Real.SingleRun_FitResults.chi2Stat.chi2min;

%% Fit
fig = figure('Renderer','opengl');
set(fig,'units','normalized','pos',[0.2, 0.2,1,0.6]);
[P,s] = polyfit(t,chi2min,1); 
ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
d=scatter(t,chi2min,'s','filled','MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));
hold on;
chi2minfit = P(1)*t+P(2);
f=plot(t,chi2minfit,'r-.','LineWidth',3,'Color',rgb('SteelBlue'));
chi2minfitu = (P(1)+ste(1))*t+P(2);
%fu=plot(t,chi2minfitu,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
chi2minfitd = (P(1)-ste(1))*t+P(2);
%fd=plot(t,chi2minfitd,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
hold off
%datetick('x',7);
xlabel('days');
xtickangle(45);
ylabel('\chi ^2 ');
title(sprintf('KATRIN - KNM1 - \\chi ^2 Stability - Start: %s - Stop: %s',...
    datestr(Real.SingleRunData.StartTimeStamp(1)),...
    datestr(Real.SingleRunData.StartTimeStamp(end)+seconds(Real.SingleRunData.TimeSec(end)))));
leg = legend([d,f],sprintf('data - %.0f runs - last -40 eV below E_0',...
    numel(t)),sprintf('linear fit - slope = %0.3f \\pm %0.3f per week',P(1)*7,ste(1)*7),...
    'Location','northeast');
leg.Color = 'none'; legend boxoff;
PrettyFigureFormat ; set(gca,'FontSize',18);
export_fig(gcf,'plots/KNM1_Chi2Stability.png','-m3');

%% Plot Chi2 distribution superimposed to a Chi2 low
 FitChi2distribution(Real.SingleRun_FitResults.chi2Stat.chi2min);
 PrettyFigureFormat ; set(gca,'FontSize',18);
export_fig(gcf,'plots/KNM1_Chi2distribution.png','-m3');

%% Runs with the lowest p-values

[Maxchi2 IndexMaxChi2] = sort(chi2min); 

% 1st
RunMaxChi2_1 = Real.RunList(IndexMaxChi2(end));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_1,Maxchi2(end));
A=RunAnalysis('RunNr',RunMaxChi2_1,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_1)),'-m3');
% 2nd
RunMaxChi2_2 = Real.RunList(IndexMaxChi2(end-1));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_2,Maxchi2(end-1));
A=RunAnalysis('RunNr',RunMaxChi2_2,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_2)),'-m3');
% 3rd
RunMaxChi2_3 = Real.RunList(IndexMaxChi2(end-2));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_3,Maxchi2(end-2));
A=RunAnalysis('RunNr',RunMaxChi2_3,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_3)),'-m3');
% 4th
RunMaxChi2_4 = Real.RunList(IndexMaxChi2(end-3));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_4,Maxchi2(end-3));
A=RunAnalysis('RunNr',RunMaxChi2_4,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_4)),'-m3');
% 5th
RunMaxChi2_5 = Real.RunList(IndexMaxChi2(end-4));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_5,Maxchi2(end-4));
A=RunAnalysis('RunNr',RunMaxChi2_5,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_5)),'-m3');
% 6th
RunMaxChi2_6 = Real.RunList(IndexMaxChi2(end-6));
fprintf('Run: %d - chi2 = %.1f\n',RunMaxChi2_6,Maxchi2(end-5));
A=RunAnalysis('RunNr',RunMaxChi2_6,'DataType','Real'); A.exclDataStart=14 ;A.Fit ; A.PlotFit('saveplot','OFF') ;
export_fig(gcf,sprintf('SingleRunFit40eVStat_Run%s.png',num2str(RunMaxChi2_6)),'-m3');




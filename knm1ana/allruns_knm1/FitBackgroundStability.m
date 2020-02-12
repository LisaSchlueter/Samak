%
% Fit Background Stability with a Constant / Constant + Slope
% T. Lasserre
% Last Modified, June 3 2019
%

% Define Data Set
t     =  (datenum(Real.SingleRunData.StartTimeStamp)-datenum(Real.SingleRunData.StartTimeStamp(1)));%days
bkg    =  Real.SingleRun_FitResults.chi2Stat.B*1000;% - mean(Real.SingleRun_FitResults.chi2Stat.B);
bkgerr =  Real.SingleRun_FitResults.chi2Stat.BErr*1000;%./mean(Real.SingleRun_FitResults.chi2Stat.B);

%% Fit
fig = figure('Renderer','opengl');
set(fig,'units','normalized','pos',[0.2, 0.2,1,0.6]);
[P,s] = polyfit(t,bkg,1); 
ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
d=scatter(t,bkg,'s','filled','MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));
hold on;
bkgfit = P(1)*t+P(2);
f=plot(t,bkgfit,'r-.','LineWidth',3,'Color',rgb('SteelBlue'));
bkgfitu = (P(1)+ste(1))*t+P(2);
%fu=plot(t,bkgfitu,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
bkgfitd = (P(1)-ste(1))*t+P(2);
%fd=plot(t,bkgfitd,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
hold off
%datetick('x',7);
xlabel('days');
xtickangle(45);
ylabel('B (mcps) ');
title(sprintf('KATRIN - KNM1 - Background Stability - Start: %s - Stop: %s',...
    datestr(Real.SingleRunData.StartTimeStamp(1)),...
    datestr(Real.SingleRunData.StartTimeStamp(end)+seconds(Real.SingleRunData.TimeSec(end)))));
leg = legend([d,f],sprintf('data - %.0f runs - last -40 eV below E_0',...
    numel(t)),sprintf('linear fit - slope = %.2f \\pm %.2f mcps per week',P(1)*7,ste(1)*7),'Location','northeast');
            leg.Color = 'none'; legend boxoff;
PrettyFigureFormat ; set(gca,'FontSize',18);
export_fig(gcf,'plots/KNM1_BackgroundStability.png','-m3');
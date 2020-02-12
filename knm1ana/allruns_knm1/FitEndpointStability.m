%
% Fit Endpoint Stability with a Constant / Constant + Slope
% T. Lasserre
% Last Modified, June 3 2019
%

% Define Data Set
t     =  (datenum(Real.SingleRunData.StartTimeStamp)-datenum(Real.SingleRunData.StartTimeStamp(1)));%days
e0    =  Real.SingleRun_FitResults.chi2Stat.E0 - mean(Real.SingleRun_FitResults.chi2Stat.E0);
e0err =  Real.SingleRun_FitResults.chi2Stat.E0Err./mean(Real.SingleRun_FitResults.chi2Stat.E0);

%% Fit
fig = figure('Renderer','opengl');
set(fig,'units','normalized','pos',[0.2, 0.2,1,0.6]);
%[P] = polyfitweighted(t,e0,1,1./e0err); 
[P,s] = polyfit(t,e0,1); 
ste = sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
d=scatter(t,e0,'s','filled','MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));
hold on;
e0fit = P(1)*t+P(2);
f=plot(t,e0fit,'r-.','LineWidth',3,'Color',rgb('SteelBlue'));
e0fitu = (P(1)+ste(1))*t+P(2);
%fu=plot(t,e0fitu,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
e0fitd = (P(1)-ste(1))*t+P(2);
%fd=plot(t,e0fitd,'r-.','LineWidth',2,'Color',rgb('CadetBlue'));
hold off
%datetick('x',7);
xlabel('days');
xtickangle(45);
ylabel('E_0 - <E_0> (eV) ');
title(sprintf('KATRIN - KNM1 - Endpoint Stability - Start: %s - Stop: %s',...
    datestr(Real.SingleRunData.StartTimeStamp(1)),...
    datestr(Real.SingleRunData.StartTimeStamp(end)+seconds(Real.SingleRunData.TimeSec(end)))));
leg = legend([d,f],sprintf('data - %.0f runs - last -40 eV below E_0',numel(t)),sprintf('linear fit - slope = %0.3f \\pm %0.3f eV per week',P(1)*7,ste(1)*7),'Location','southeast');
            leg.Color = 'none'; legend boxoff;
PrettyFigureFormat ; set(gca,'FontSize',18);
export_fig(gcf,'plots/KNM1_EndpointStability.png','-m3');
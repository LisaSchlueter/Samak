function m2sysVsyear2019()
% Plot of Squared neutrino mass stat/syst errors obtained from tritium beta decay
% in the  period 1990-2019 plotted against the year of publication 
% Previous: results from the more recent experiments in Mainz and Troitsk
% New: results from KATRIN KNM1

%% Load Data
[mBetaSquared , mBetaSquaredPDG]  = mBetaSquaredExperiments();

%% Plot
maintitle=sprintf('Squared neutrino mass Uncertainties obtained from tritium \\beta -decay in the period 1990-2019');
savefile=sprintf('plots/m2sysVsyear2019.png');
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1400 600]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=26;a.FontWeight='normal';
stem([mBetaSquared.Year], [mBetaSquared.Sys],'LineStyle','--','Color','Black','LineWidth',2,'Marker','none');
hold on
for i=1:12
    h(i)=errorbar(mBetaSquared(i).Year, mBetaSquared(i).Sys,mBetaSquared(i).Sys*0,'Marker',mBetaSquared(i).Marker,'MarkerSize',12,'LineWidth',1,'Color','Black',...
     'MarkerEdgeColor',mBetaSquared(i).Color,'MarkerFaceColor',mBetaSquared(i).Color);
     h(i).CapSize = 0; h(i).LineStyle= 'none';h(i).LineWidth= 3;
    hold on
end
l=legend(h,[mBetaSquared(1).Experiment ' ' mBetaSquared(2).Reference],...
    [mBetaSquared(2).Experiment ' ' mBetaSquared(2).Reference],...
       [mBetaSquared(3).Experiment ' ' mBetaSquared(3).Reference],...
    [mBetaSquared(4).Experiment ' ' mBetaSquared(4).Reference],...
    [mBetaSquared(5).Experiment ' ' mBetaSquared(5).Reference],...
    [mBetaSquared(6).Experiment ' ' mBetaSquared(6).Reference],...
    [mBetaSquared(7).Experiment ' ' mBetaSquared(7).Reference],...
    [mBetaSquared(8).Experiment ' ' mBetaSquared(8).Reference],...
    [mBetaSquared(9).Experiment ' ' mBetaSquared(9).Reference],...
    [mBetaSquared(10).Experiment ' ' mBetaSquared(10).Reference],...
    [mBetaSquared(11).Experiment ' ' mBetaSquared(11).Reference],...
    [mBetaSquared(12).Experiment ' ' mBetaSquared(12).Reference],...
    'FontSize',12,'Location','north');
l.NumColumns=4;
legend('boxoff');
xlabel('year');
ylabel('m^2(\nu_e) c^4 systematic uncertainty (eV^2)');
grid on
FontName = 'Arial';
FontSize = 22;
set(gca,'FontName',FontName,'FontSize',FontSize);
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize);
set(get(gca,'YLabel'),'FontName',FontName,'FontSize',FontSize);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','on');
set(gcf,'Color','w');
set(gca,'YScale','log');
box on;
xlim([mBetaSquared(1).Year-1 mBetaSquared(12).Year+1]);
ylim([1e-1 3*max([mBetaSquared.Sys])])
export_fig(gcf,savefile,'-r300');
%
% Given a Run Number
% Fit of the data & MC
% Results Superimposed
% and compared
% 
% Th. Lasserre
% 20/05/2019
% 

% Inputs
Run           = 51458;
exclDataStart = 14;
qUPlotShift   = 0.5; 

% Read Data / Twin-Asimov
R=RunAnalysis('RunNr',Run,'DataType','Real','exclDataStart',exclDataStart);
T=RunAnalysis('RunNr',Run,'DataType','Twin','exclDataStart',exclDataStart);

% Fit Real / Twin-Asimov
R.Fit; 
T.Fit;

%% Plot Results
R.PlotFit;
T.PlotFit

%% SuperImpose Results
fig9 = figure(9);
% Spectrum
s1= subplot(2,1,1);
mR = plot(R.ModelObj.qU-R.ModelObj.Q_i,R.ModelObj.TBDIS./R.ModelObj.qUfrac./R.ModelObj.TimeSec,...
    'Color',rgb('Amethyst'));
hold on;
dR = errorbar(R.RunData.qU-R.ModelObj.Q_i,R.RunData.TBDIS./R.ModelObj.qUfrac./R.ModelObj.TimeSec,R.RunData.TBDISE./R.ModelObj.qUfrac./R.ModelObj.TimeSec,...
    'd','MarkerSize',5,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2,'Color',rgb('Amethyst'));
dT = errorbar(T.RunData.qU-T.ModelObj.Q_i,T.RunData.TBDIS./T.ModelObj.qUfrac./T.ModelObj.TimeSec,T.RunData.TBDISE./T.ModelObj.qUfrac./T.ModelObj.TimeSec,...
    'o','MarkerSize',5,'MarkerFaceColor',rgb('Orange'),'LineWidth',2,'Color',rgb('Orange'));
hold off;
myleg = legend([mR dR dT],'Data Best Fit','Data','Twin MC','Location','Northeast');
legend('boxoff');
xlabel(sprintf('retarding potential - %.1f (eV)',R.ModelObj.Q_i),'FontSize',16);
ylabel('cps');
set(gca,'yscale','log');
set(gca,'YMinorTick','off');
set(gca,'TickLength',[0.01 0.01]);
xlim([min(min(R.ModelObj.qU(R.exclDataStart)))-R.ModelObj.Q_i max(max(R.ModelObj.qU))-R.ModelObj.Q_i]);
ylim([0 1.25*max(R.ModelObj.TBDIS(R.exclDataStart)./R.ModelObj.qUfrac(R.exclDataStart)./R.ModelObj.TimeSec)]);
PrettyFigureFormat;
                title(sprintf('Run%u: Data and Twin MonteCarlo (Asimov)',Run));
set(gca,'FontSize',22);

% Residual wrt Best Fit
s2= subplot(2,1,2);
DataMinusBestFitE = R.RunData.TBDISE.*0;
DataMinusBestFit  = (R.RunData.TBDIS-R.ModelObj.TBDIS)./R.RunData.TBDISE;
TwinMinusBestFitE = T.RunData.TBDISE.*0;
TwinMinusBestFit  = (T.RunData.TBDIS-R.ModelObj.TBDIS)./T.RunData.TBDISE;
dR = errorbar(R.RunData.qU-R.ModelObj.Q_i-qUPlotShift,DataMinusBestFit,DataMinusBestFitE,...
    'd','MarkerSize',8,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2,'Color',rgb('Amethyst'));
hold on
dT = errorbar(T.RunData.qU-T.ModelObj.Q_i+qUPlotShift,TwinMinusBestFit,TwinMinusBestFitE,...
    'o','MarkerSize',8,'MarkerFaceColor',rgb('Orange'),'LineWidth',2,'Color',rgb('Orange'));
hold off
xlim([min(min(R.ModelObj.qU(R.exclDataStart)))-R.ModelObj.Q_i max(max(R.ModelObj.qU))-R.ModelObj.Q_i]);
myleg = legend([dR dT],'Data - Best Fit Data','Twin MC - Best Fit Data','Location','Northeast');
legend('boxoff');
PrettyFigureFormat;
xlabel(sprintf('retarding potential - %.1f (eV)',R.ModelObj.Q_i),'FontSize',16);
ylabel('Normalized Residuals');
set(gca,'FontSize',22);

savefile = ['./plots/' 'RealTwin_RunFit' num2str(Run) '.png'];
export_fig(fig9,savefile);

%% settings
RunList               = 'KNM1';
exclDataStart         = 2;
nTrials               = 5000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
BkgCM                 = 'ON';

%% Init Model Object and covariance matrix object
TwinALLON = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.1);
TwinALLON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');

%% Sensitivity - Via ASimov Fit Error
%for index=2:25
%%
counter=0;
for index=2:25
    counter=counter+1;
    TwinALLON.exclDataStart=index;
    TwinALLON.NonPoissonScaleFactor=1.1;
    % Stat
    TwinALLON.chi2='chi2Stat';TwinALLON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
    TwinALLON.Fit; %TwinALLON.PlotFit;
    Rstat = TwinALLON.FitResult;
    Smsquared(counter)     = Rstat.par(1);
    Smsquared_err(counter) = Rstat.err(1);
    % CM Shape
%     TwinALLON.NonPoissonScaleFactor=1.1;
%     TwinALLON.chi2='chi2CMShape';TwinALLON.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM); 
%     TwinALLON.Fit; %TwinALLON.PlotFit;
%     Rsysshape = TwinALLON.FitResult;
%     SSmsquared(counter)     = Rsysshape.par(1);
%     SSmsquared_err(counter) = Rsysshape.err(1);
end

%% KNM1 Stat / Syst
T     = days(seconds(cumsum(TwinALLON.RunData.TimeSec)));
index=2:25;
counter=1:24;
maintitle=sprintf('KNM1 Neutrino Mass Square Uncertainty - %.1f days - Golden Run List - Golden Pixel List',T);
savefile=sprintf('plots/KNM1NeutrinoMassSquareUncertainty%.0fdaysGRPL-1.png',T);
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

h2=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i, SSmsquared_err(counter),'LineWidth',7);
hold on
h1=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i,Smsquared_err(counter),'LineWidth',7)
hold off
leg = legend([h1,h2],'stat.','stat. + syst. (shape-only)','Location','northwest');
leg.Color = 'none'; legend boxoff;
xlabel(sprintf('Retarding Potential - %.1f (eV)',TwinALLON.ModelObj.Q_i));
ylabel('Neutrino Mass Squared Error (eV^2)');
grid on
PrettyFigureFormat
set(gca,'FontSize',24);
export_fig(gcf,savefile,'-q101','-m3');

%% KNM1 Stat / Syst Breakdown
maintitle=sprintf('KNM1 Neutrino Mass Square Uncertainty - %.1f days - Golden Run List - Golden Pixel List',T);
savefile=sprintf('plots/KNM1NeutrinoMassSquareUncertainty%.0fdaysGRPL-2.png',T);
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
h2=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i, sqrt(SSmsquared_err(counter).^2-Smsquared_err(counter).^2),'LineWidth',7);
hold on
h1=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i,Smsquared_err(counter),'LineWidth',7);
%h3=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i, (SSmsquared_err(index)),'LineWidth',5,'LineStyle','-','Color','Black');
hold off
leg = legend([h1,h2],'stat.','syst. (shape-only)','Location','northwest');
leg.Color = 'none'; legend boxoff;
xlabel(sprintf('Retarding Potential - %.1f (eV)',TwinALLON.ModelObj.Q_i));
ylabel('Neutrino Mass Squared Error (eV^2)');
grid on
PrettyFigureFormat
set(gca,'FontSize',24);
export_fig(gcf,savefile,'-q101','-m3');

%% KNM1 Asimov Sensitivity + Stat / Syst Breakdown
maintitle=sprintf('KNM1 Neutrino Mass Square Uncertainty - %.1f days - Golden Run List - Golden Pixel List',T);
savefile=sprintf('plots/KNM1NeutrinoMassSquareUncertainty%.0fdaysGRPL-3.png',T);
fig1 = figure('Name','Test','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
h1=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i, SSmsquared_err(counter),'LineWidth',10);
hold on
h2=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i,Smsquared_err(counter),'LineWidth',5);
h3=stairs(TwinALLON.RunData.qU(index)-TwinALLON.ModelObj.Q_i, sqrt(SSmsquared_err(counter).^2-Smsquared_err(counter).^2),'LineWidth',5);
hold off
leg = legend([h1,h2, h3],'stat+syst','stat','syst = [(stat+syst)^2-stat^2]^{0.5}','Location','northwest');
leg.Color = 'none'; legend boxoff;
xlabel(sprintf('Retarding Energy - %.1f (eV)',TwinALLON.ModelObj.Q_i));
ylabel('Neutrino Mass Squared Error (eV^2)');
grid on
xlim([-93 -20])
PrettyFigureFormat
set(gca,'FontSize',24);
export_fig(gcf,savefile,'-q101','-m3');

%% Latex Display
t = PrintTable(sprintf('Fit Results  - Fitter = %s',TwinALLON.fitter));
t.addRow('Fit Range','$m^2$ shift','stat+syst','stat','syst','Unit');
counter=0;
for index=2:25
    counter=counter+1;
t.addRow(sprintf('$[%.1f - %.1f]$',TwinALLON.RunData.qU(index),TwinALLON.RunData.qU(end)),...
    sprintf('%.2f',SSmsquared(counter)),...
    sprintf('%.2f',SSmsquared_err(counter)),...
    sprintf('%.2f',Smsquared_err(counter)),...
    sprintf('%.2f',sqrt(SSmsquared_err(counter)^2-Smsquared_err(counter)^2)),...
    'eV$^2$');
end
t.display;
t.HasHeader = true;
t.Format = 'tex';
t.Caption = sprintf('Fit Results - Fitter = %s - 1 $\\sigma$ uncertainties',TwinALLON.fitter);
t.print;
                
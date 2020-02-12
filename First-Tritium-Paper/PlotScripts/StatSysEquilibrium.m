% FT paper plot
% display endpoint sensitivity as a function of fit range
% find equilibrium between stat and sys for optimal fit range
% Lisa August 2019

%% Set up model
RunList = 'FTpaper';
FT = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataEffCor','RunSummary',...
    'ELossFlag','Abdurashitov','exclDataStart',12,'SysBudget',0,'DataType','Twin','FSDFlag','SAENZ',...
    'fixPar','1 5 6 7 8 9 10 11','TwinBias_Q',18574.5);
FT.NonPoissonScaleFactor = 1;
qUmin = 8:1:18;

savedir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
savefile = [savedir,sprintf('FT_FitRangeStatSys_%s.mat',FT.DataType)];%&_OnlyELoss

if exist(savefile,'file')
    load(savefile);
else
%% Calculate Sensitivites

% stat
FT.chi2 = 'chi2Stat';
E0Stat = zeros(numel(qUmin),1);
progressbar('stat');
for i=1:numel(qUmin)
    progressbar(i/numel(qUmin));
    FT.exclDataStart = qUmin(i);
    FT.Fit;
    E0Stat(i) = FT.FitResult.err(2);
end

defaultEffects = struct(...
    'RF_EL','OFF',...   % Response Function(RF) EnergyLoss
    'RF_BF','ON',...   % RF B-Fields
    'RF_RX','ON',...   % Column Density, inel cross ection
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON',...
    'DOPoff','OFF',...
    'Stack','ON',...
    'FPDeff','ON');

% stat + syst
FT.chi2 = 'chi2CMShape';
progressbar('stat + sys');
E0Tot = zeros(numel(qUmin),1);
for i=1:numel(qUmin)
     progressbar(i/numel(qUmin));
    FT.ComputeCM('SysEffects',defaultEffects);
    FT.exclDataStart = qUmin(i);
    FT.Fit;
    E0Tot(i) = FT.FitResult.err(2);
end

% sys
E0Sys = (E0Tot.^2-E0Stat.^2).^0.5;
save(savefile,'E0Stat','E0Tot','E0Sys','qUmin');
end
%% plot
f11 = figure('Renderer','painters');
set(f11, 'Units', 'centimeters', 'Position', [0.1, 0.1, 8.4, 5]); % 0.8,0.75

qU = FT.RunData.qU(qUmin)-18573.7;
qUEqual = -96.5;

PlotArg = {'LineWidth',1.5,'Color'};
plotqU = linspace(min(qU),max(qU),100);
peq  = plot(qUEqual.*[1,1],[0+0.053,max(E0Tot)-0.001],'-','LineWidth',7,'Color',rgb('LightGray'));
hold on;
%pstat1 = plot(plotqU,interp1(qU,E0Stat,plotqU,'spline'),PlotArg{:},rgb('DodgerBlue'),'LineStyle','-.');
pstat = plot(plotqU,smooth(interp1(qU,E0Stat,plotqU,'lin')),PlotArg{:},rgb('DodgerBlue'),'LineStyle','--','LineWidth',1.5);

ptot = plot(plotqU,smooth(interp1(qU,E0Tot,plotqU,'lin'),10),PlotArg{:},rgb('FireBrick'),'LineStyle','-');
psys = plot(plotqU,smooth(interp1(qU,E0Sys,plotqU,'lin'),10),PlotArg{:},rgb('Orange'),'LineStyle',':');
xlabel(sprintf('Lower fit boundary below {\\itE}_0 (eV)'))
ylabel({sprintf('1\\sigma uncertainty');'of fitted endpoint (eV)'});  %1-sigma effective endpoint uncertainty 
leg = legend([pstat,psys,ptot,peq],' Statistical',' Systematic',' Total',' Statictical = systematic');
legend boxoff; leg.Location = 'northwest';
FTpaperFormat;
LocalFontSize =7;% 37.5;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+2);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+2);
leg.FontSize = LocalFontSize+2;
xlim([min(qU),-50]);%max(qU)]);
ylim([min(E0Stat),max(E0Tot)]);
set(gca,'LineWidt',0.7);
%% remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)+0.01;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4)- ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% save plot
plotfile = [strrep(strrep(savefile,'results','plots'),'.mat','.pdf')];
export_fig(f11,plotfile)
%print(f11,plotfile,'-r450','-dpng')





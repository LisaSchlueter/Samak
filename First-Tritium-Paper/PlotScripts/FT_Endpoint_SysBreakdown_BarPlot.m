% Quick Neutrino Mass Sensitivity Study for First Tritium
% Compute and Display
RunList = 'FTpaper';
RecomputeFlag = 'OFF';
Mode = 'Data';
exclDataStart = 13;
switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7 
        belowE0 = 400;
    case 9
        belowE0 = 200;
    case 13
        belowE0 = 100;
    case 12
        belowE0 = 125;
end
M = MultiRunAnalysis('RunList',RunList,'exclDataStart',exclDataStart,'DataType','Real',...
    'FSDFlag','SAENZ','ELossFlag','Abdurashitov','SysBudget',0,'DataEffCor','RunSummary',...
    'NonPoissonScaleFactor',1,'fixPar','1 5 6 7 8 9 10 11');
M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
TBDIS = M.ModelObj.TBDIS;
TBDISE = M.ModelObj.TBDISE;
chi2name = 'chi2CMShape';
savedir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
savename = [savedir,sprintf('FT_EndpointSysBreakdownData_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
d = importdata(savename);
E0StackedData = d.E0Stacked;
E0SingleData  = d.E0Single;
E0MultiData = d.E0Multi;
M.chi2 = 'chi2Stat';
savename = [savedir,sprintf('FT_EndpointSensitivity_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
load(savename);

%% get NP component
% Data Stack
M.NonPoissonScaleFactor = 1.03;
M.chi2 = 'chi2CMShape';
M.ComputeCM;
M.Fit;
E0StackedData = [E0StackedData;M.FitResult.err(2)-E0MultiData(end)];
% Data Single
M.chi2 = 'chi2Stat';
M.Fit;
E0SingleData = [E0SingleData;M.FitResult.err(2)-E0SingleData(1)];
%

% Sim 
M.RunData.qU = M.ModelObj.qU;
M.RunData.TBDIS = TBDIS;
M.RunData.TBDISE = TBDISE;
M.NonPoissonScaleFactor = 1.03;

% sim stack
M.chi2 = 'chi2CMShape';
M.ComputeCM;
M.Fit;
E0Stacked = [E0Stacked;M.FitResult.err(2)-E0Multi(end)];
% sim single
M.chi2 = 'chi2Stat';
M.Fit;
E0Single = [E0Single;M.FitResult.err(2)-E0Single(1)];
%%

%% Display
f55 = figure('Name','MultiBarPlot','Renderer','painters');
set(f55, 'Units', 'centimeters', 'Position', [0.1, 0.1, 8.4, 12]);
b = barh([10;20],[E0Stacked , E0StackedData]','stacked');
b(1).LineStyle = 'none';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';
b(6).LineStyle = 'none';b(7).LineStyle = 'none'; b(8).LineStyle = 'none';b(9).LineStyle = 'none'; b(10).LineStyle = 'none';
b(1).FaceColor = rgb('Gainsboro');
b(3).FaceColor = rgb('FireBrick');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('Navy');
b(7).FaceColor = rgb('PowderBlue');
b(8).FaceColor = rgb('YellowGreen');
b(9).FaceColor = rgb('DarkGreen');
b(10).FaceColor = rgb('DarkGreen');
leg_str = {'Statistical uncertainty';'+ DT concentration fluctuation';'+ Final-state distribution';...
    '+ Magnetic fields';'+ Energy-loss function'; sprintf('+ Column density, Inel. cross-section');...
    '+ Detector efficiency';'+ Background slope';'Non-Poisson background'}; %;'+ Response Function'};
b(1).BarWidth = b(1).BarWidth/2;
%% Single Sys
SingleSysBarWidth = b(1).BarWidth/(2*3);
hold on;
pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
%bTC  = barh([7.7;21], [mNuSingle(1),mNuSingle(2); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
%bTC(2).FaceColor =rgb('FireBrick'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
%bTC(1).FaceColor = 'w';%bFSD(1).FaceColor = 'w'; bRF(1).FaceColor = 'w';

bTASR  = barh([7.65;17.7], [E0Single(1),E0Single(3); E0SingleData(1),E0SingleData(3)],'stacked','BarWidth',SingleSysBarWidth);
bTASR(2).FaceColor =rgb('FireBrick'); bTASR(2).LineStyle = 'none'; bTASR(2).FaceAlpha = 0.5;
bTASR(1).FaceColor = 'w'; bTASR(1).LineStyle = 'none';

bFSD  = barh([7;17.05], [E0Single(1),E0Single(4); E0SingleData(1),E0SingleData(4)],'stacked','BarWidth',SingleSysBarWidth);
bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
bFSD(1).FaceColor = 'w'; bFSD(1).LineStyle = 'none';

bRF_BF  = barh([6.35;16.4], [E0Single(1),E0Single(5); E0SingleData(1),E0SingleData(5)],'stacked','BarWidth',SingleSysBarWidth);
bRF_BF(2).FaceColor =rgb('CadetBlue'); bRF_BF(2).LineStyle = 'none'; bRF_BF(2).FaceAlpha = 0.5;
bRF_BF(1).FaceColor = 'w'; bRF_BF(1).LineStyle = 'none';

bRF_EL  = barh([5.7;15.75], [E0Single(1),E0Single(6);  E0SingleData(1),E0SingleData(6)],'stacked','BarWidth',SingleSysBarWidth);
bRF_EL(2).FaceColor =rgb('PowderBlue'); bRF_EL(2).LineStyle = 'none'; bRF_EL(2).FaceAlpha = 0.9;
bRF_EL(1).FaceColor = 'w'; bRF_EL(1).LineStyle = 'none';

bRF_RX  = barh([5.05;15.1], [E0Single(1),E0Single(7); E0SingleData(1),E0SingleData(7)],'stacked','BarWidth',SingleSysBarWidth);
bRF_RX(2).FaceColor =rgb('DarkSlateGray'); bRF_RX(2).LineStyle = 'none'; bRF_RX(2).FaceAlpha = 0.5;
bRF_RX(1).FaceColor = 'w'; bRF_RX(1).LineStyle = 'none';

bFPD  = barh([4.4;14.44], [E0Single(1),E0Single(8); E0SingleData(1),E0SingleData(8)],'stacked','BarWidth',SingleSysBarWidth);
bFPD(2).FaceColor =rgb('YellowGreen'); bFPD(2).LineStyle = 'none'; bFPD(2).FaceAlpha = 0.5;
bFPD(1).FaceColor = 'w'; bFPD(1).LineStyle = 'none';

bBkg  = barh([3.74;13.77], [E0Single(1),E0Single(9); E0SingleData(1),E0SingleData(9)],'stacked','BarWidth',SingleSysBarWidth);
bBkg(2).FaceColor =rgb('DarkGreen'); bBkg(2).LineStyle = 'none'; bBkg(2).FaceAlpha = 0.5;
bBkg(1).FaceColor = 'w'; bBkg(1).LineStyle = 'none';

bBkgNP  = barh([3.09;13.1], [E0Single(1),E0Single(10); E0SingleData(1),E0SingleData(10)],'stacked','BarWidth',SingleSysBarWidth);
bBkgNP(2).FaceColor =rgb('HotPink'); bBkgNP(2).LineStyle = 'none'; bBkgNP(2).FaceAlpha = 0.5;
bBkgNP(1).FaceColor = 'w'; bBkgNP(1).LineStyle = 'none';
%% Plot again for nicer statistical border --
b = barh([10;20],[E0Stacked , E0StackedData]','stacked');
b(1).LineStyle = ':'; b(1).LineWidth = 1.0;b(1).EdgeColor = rgb('DimGray');
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';
b(5).LineStyle = 'none'; b(6).LineStyle = 'none';b(7).LineStyle = 'none';
b(8).LineStyle = 'none';b(9).LineStyle = 'none';b(10).LineStyle = 'none';
b(1).FaceColor = rgb('Gainsboro');
b(3).FaceColor = rgb('FireBrick');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('PowderBlue');
b(7).FaceColor = rgb('DarkSlateGray');
b(8).FaceColor = rgb('YellowGreen');
b(9).FaceColor = rgb('DarkGreen');
b(10).FaceColor = rgb('DarkGreen');
b(1).BarWidth = b(1).BarWidth/2;

leg_str = {leg_str{:},'Single Contributions','Stat + Tritium Acitivity Fluctuation', 'Stat + Final State Distribution',...
    'Stat + Magnetic Fields','Stat + Energy Loss Function','Stat + Column Density + Inel. Cross Section',...
    'Stat + FPD Efficiency','Stat + Background Slope'};
leg = legend([b(1),b(3),b(4),b(5),b(6),b(7),b(8),b(9),pnone,bTASR(2),bFSD(2),bRF_BF(2),bRF_EL(2),bRF_RX(2),bFPD(2),bBkg(2)],leg_str{:});
leg.Location = 'north';
leg.NumColumns = 2; legend boxoff

FTpaperFormat;
ylim([2.2 33]);
xlim([0.15,0.27]);
set(gca,'FontSize',9);

Pos_tmp = get(gcf,'Position');

if Pos_tmp(3)<10
    leg = legend([b(1),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10)],cellfun(@(x) strrep(x,'+ ',''),leg_str(1:9),'UniformOutput',0));
    leg.NumColumns=1;
    leg.FontSize=8;
    Pos2_tmp = leg.Position;
    leg.Position= [Pos2_tmp(1)+0.025,Pos2_tmp(2)+0.01,Pos2_tmp(3:4)];
else
    leg.FontSize = 7;
    mypos = leg.Position;
    leg.Position = [mypos(1)+0.017,mypos(2),mypos(3:4)];
end

xlabel(sprintf('1\\sigma uncertainty of fitted endpoint (eV)'));
yticks = [10,20];yticklabels({'Simulation','Data'});
xticks(0:0.02:0.28)
set(gca,'YMinorTick','off');
h = gca;
 h.YRuler.TickLength = [0,0];
%set(gca,'TickLength',[.01 0]);
%%
% leg.delete
% set(f55, 'Units', 'centimeters', 'Position', [0.1, 0.1, 8.4, 5]);
% ylim([2 24]);
%% Save
plotdir =strrep(savedir,'results','plots');
MakeDir(plotdir);
saveplot = [plotdir,sprintf('FT_E0SysBreakdown_%.0feV_%s_%s.pdf',belowE0,chi2name,RunList)];
export_fig(gcf,saveplot);








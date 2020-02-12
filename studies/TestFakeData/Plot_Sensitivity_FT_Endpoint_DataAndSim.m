% Quick Neutrino Mass Sensitivity Study for First Tritium
% Compute and Display
%RecomputeFlag = 'OFF';
%Mode = 'Data';
exclDataStart = 13;
RunList = 'FTpaper';
switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7
        belowE0 = 400;
    case 9
        belowE0 = 200;
    case 12
        belowE0 = 125;
    case 13
        belowE0 = 100;
end

chi2name = 'chi2CMShape';
%mydir = [getenv('SamakPath'),'studies/local_LisaThesisPlots/results/'];
mydir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
savename = [mydir,sprintf('FT_EndpointSysBreakdownData_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
d = importdata(savename);
E0StackedData = d.E0Stacked;
E0SingleData  = d.E0Single;

savename =[mydir,sprintf('FT_EndpointSensitivity_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
load(savename);

%% Display
f55 = figure('Name','MultiBarPlot','Renderer','painters');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
b = barh([10;20],[E0Stacked , E0StackedData]','stacked');
b(1).LineStyle = 'none';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';
b(6).LineStyle = 'none';b(7).LineStyle = 'none';b(8).LineStyle = 'none';b(9).LineStyle = 'none';
leg_str = {' Statistical';'+ Tritium Acitivity Fluctuation';'+ Final State Distribution';...
    '+ Magnetic Fields';'+ Energy Loss Function'; '+ Column Density + Inel. Cross Section';'+ FPD Efficiency';'+ Background Slope'}; %;'+ Response Function'};
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

bRF_RX  = barh([5.03;15.08], [E0Single(1),E0Single(7); E0SingleData(1),E0SingleData(7)],'stacked','BarWidth',SingleSysBarWidth);
bRF_RX(2).FaceColor =rgb('DarkSlateGray'); bRF_RX(2).LineStyle = 'none'; bRF_RX(2).FaceAlpha = 0.5;
bRF_RX(1).FaceColor = 'w'; bRF_RX(1).LineStyle = 'none';

bFPD  = barh([4.36;14.43], [E0Single(1),E0Single(8); E0SingleData(1),E0SingleData(8)],'stacked','BarWidth',SingleSysBarWidth);
bFPD(2).FaceColor =rgb('YellowGreen'); bFPD(2).LineStyle = 'none'; bFPD(2).FaceAlpha = 0.5;
bFPD(1).FaceColor = 'w'; bFPD(1).LineStyle = 'none';

bBkg  = barh([3.69;13.76], [E0Single(1),E0Single(9); E0SingleData(1),E0SingleData(9)],'stacked','BarWidth',SingleSysBarWidth);
bBkg(2).FaceColor =rgb('DarkGreen'); bBkg(2).LineStyle = 'none'; bBkg(2).FaceAlpha = 0.5;
bBkg(1).FaceColor = 'w'; bBkg(1).LineStyle = 'none';
%% Plot again for nicer statistical border --
b = barh([10;20],[E0Stacked , E0StackedData]','stacked');
b(1).LineStyle = '--'; b(1).LineWidth=2.5;
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none'; b(6).LineStyle = 'none';
b(7).LineStyle = 'none';b(8).LineStyle = 'none';b(9).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(3).FaceColor = rgb('FireBrick');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('PowderBlue');
b(7).FaceColor = rgb('DarkSlateGray');
b(8).FaceColor = rgb('YellowGreen');
b(9).FaceColor = rgb('DarkGreen');
b(1).BarWidth = b(1).BarWidth/2;

leg_str = {leg_str{:},'Single Contributions','Stat + Tritium Acitivity Fluctuation', 'Stat + Final State Distribution',...
    'Stat + Magnetic Fields','Stat + Energy Loss Function','Stat + Column Density + Inel. Cross Section','Stat + FPD Efficiency','Stat + Background Slope'};
leg = legend([b(1),b(3),b(4),b(5),b(6),b(7),b(8),b(9),pnone,bTASR(2),bFSD(2),bRF_BF(2),bRF_EL(2),bRF_RX(2),bFPD(2),bBkg(2)],leg_str{:});
leg.Location = 'north';
leg.NumColumns = 2; legend boxoff
FTpaperFormat;
ylim([2.5 33]);
xlim([0.15,0.27]);
set(gca,'FontSize',29);
leg.FontSize=22;
%xlabel(sprintf('endpoint uncertainty / sensitivity 68.3%% C.L. (eV)'));
xlabel(sprintf('1\\sigma effective endpoint uncertainty (eV)'));
yticks = [10,20];yticklabels({'Simulation','Data'});
xticks(0:0.02:0.28)
 %% remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)+0.01;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.012;
ax_height = outerpos(4)- ti(2) - ti(4)-0.001;
ax.Position = [left bottom ax_width ax_height];
%% Save
plotdir = strrep(mydir,'results','plots');
plot_name = [plotdir,sprintf('FT_E0SensitivityBreakdown_Comparison_%.0feVbelowE0.png',belowE0)];
%print(f55,plot_name,'-dpng');
publish_figurePDF(f55,strrep(plot_name,'.png','.pdf'));







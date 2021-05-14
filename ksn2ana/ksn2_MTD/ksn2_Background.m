% plot different contours based on simulation with differnt background
% other settings: KNM-2 like
% look at sensitivity contour from KNM-2 like simulation with different MTDs

FakeInitFile = @ref_KNM2_KATRIN_RegMTD;
range = 40;
freePar = 'E0 Norm Bkg';
nGridSteps = 25;

%% tritium run model
F = RunAnalysis('RunNr',1,...
    'DataType','Fake',...
    'FakeInitFile',FakeInitFile,...
    'fixPar',freePar,...
    'SysBudget',40,...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'DopplerEffectFlag','FSD',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos');

F.exclDataStart = F.GetexclDataStart(range);

%% configure Sterile analysis object
SterileArg = {'RunAnaObj',F,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
%% regular KNM-2 MTD
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_RegMTD;
S.LoadGridFile
S.Interp1Grid('Maxm4Sq',34^2);
preg = S.ContourPlot('Color',rgb('DodgerBlue'),'LineStyle','-','HoldOn','OFF');
mNu4Sq_contour_reg = S.mNu4Sq_contour;
sin2T4_contour_reg = S.sin2T4_contour;

% 0 mcps
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_RegMTD_Bkg0mcps;
S.LoadGridFile
S.Interp1Grid('Maxm4Sq',34^2);
p0 = S.ContourPlot('Color',rgb('FireBrick'),'LineStyle','-.','HoldOn','ON');
mNu4Sq_contour_0 = S.mNu4Sq_contour;
sin2T4_contour_0 = S.sin2T4_contour;

% 10 mcps
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_RegMTD_Bkg10mcps;
S.LoadGridFile
S.Interp1Grid('Maxm4Sq',34^2);
p10 = S.ContourPlot('Color',rgb('ForestGreen'),'LineStyle',':','HoldOn','ON');
mNu4Sq_contour_10 = S.mNu4Sq_contour;
sin2T4_contour_10 = S.sin2T4_contour;

% 100 mcps
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_RegMTD_Bkg100mcps;
S.LoadGridFile
S.Interp1Grid('Maxm4Sq',34^2);
p100 = S.ContourPlot('Color',rgb('Orange'),'LineStyle','--','HoldOn','ON');
mNu4Sq_contour_109 = S.mNu4Sq_contour;
sin2T4_contour_109 = S.sin2T4_contour;
% change title and legend
if ~contains(freePar,'mNu')
    title(sprintf('Simulation , {\\itm}_\\nu^2 = 0 eV^2'),'FontSize',get(gca,'FontSize'),'FontWeight','normal');
else
    title(sprintf('Simulation , {\\itm}_\\nu^2 free'),'FontSize',get(gca,'FontSize'),'FontWeight','normal');
end
leg = legend([preg,p100,p10,p0],'220 mcps (KNM-2)','100 mcps','10 mcps','0 mcps');
xlim([5e-03,0.5])
leg.Title.String = sprintf('Background rate');
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);

%% save plot
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_MTD/plots/'];
MakeDir(plotdir);
plotname = sprintf('%sksn2_SensitivityCompareMTD_%s.png',plotdir,strrep(freePar,' ',''));
print(plotname,'-dpng','-r300');
fprintf('save plot to %s\n',plotname);

%%

mNu4Sq = linspace(max([min(mNu4Sq_contour_10),min(mNu4Sq_contour_reg)]),min([max(mNu4Sq_contour_10),max(mNu4Sq_contour_reg)]),1e3);
sin2T4_reg = interp1(mNu4Sq_contour_reg,sin2T4_contour_reg,mNu4Sq,'spline');
sin2T4_10 = interp1(mNu4Sq_contour_10,sin2T4_contour_10,mNu4Sq,'spline');
f2 = figure('Units','normalized','Position',[-0.1,0.1,0.6,0.6]);
s1 = subplot(2,3,[1,2,4,5]);
preg = plot(sin2T4_reg,mNu4Sq,'LineWidth',2.5);
hold on;
p10 = plot(sin2T4_10,mNu4Sq,'-.','LineWidth',2.5);
set(gca,'YScale','log')
set(gca,'XScale','log')
PrettyFigureFormat('FontSize',22);
ax1 = gca;
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV)'));
leg = legend([preg,p10],'220 mcps (KNM-2)','10 mcps');
PrettyLegendFormat(leg);
xlim([5e-03 0.5]);

s2= subplot(2,3,[3,6]);
plot(sin2T4_10.^2./sin2T4_reg.^2,mNu4Sq,'-k','LineWidth',2.5)
set(gca,'YScale','log')
PrettyFigureFormat('FontSize',22);
ax2 = gca;
ax2.YAxisLocation = 'right';
xlabel(sprintf('\\sigma_{10 mcps}^2 / \\sigma_{220 mcps}^2'));
xlim([0 1]);

linkaxes([s1,s2],'y');
ylim([1.25 max(ylim)]);
ax1.Position(3) = 0.55;
ax1.Position(2) = 0.15;
ax2.Position(2) = ax1.Position(2);
ax1.Position(1) = 0.1;
ax2.Position(1) = 0.66;
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
S.Interp1Grid;
preg = S.ContourPlot('Color',rgb('DodgerBlue'),'LineStyle','-','HoldOn','OFF');
% flat KNM-2 MTD
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_FlatMTD;
S.LoadGridFile
S.Interp1Grid;
pflat = S.ContourPlot('Color',rgb('Orange'),'LineStyle','-.','HoldOn','ON');
% isostat KNM-2 MTD
S.RunAnaObj.FakeInitFile = @ref_KNM2_KATRIN_IsoStatMTD;
S.LoadGridFile
S.Interp1Grid;
piso = S.ContourPlot('Color',rgb('ForestGreen'),'LineStyle',':','HoldOn','ON');
%% change title and legend
if ~contains(freePar,'mNu')
    title(sprintf('Simulation , {\\itm}_\\nu^2 = 0 eV^2'),'FontSize',get(gca,'FontSize'),'FontWeight','normal');
else
    title('Simulation','FontSize',get(gca,'FontSize'),'FontWeight','normal');
end

leg = legend([preg,piso,pflat],'Regular KNM-2','Iso-Stat','Flat');
xlim([5e-03,0.5])
leg.Title.String = sprintf('MTD');
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);

%% save plot
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_MTD/plots/'];
MakeDir(plotdir);
plotname = sprintf('%sksn2_SensitivityCompareMTD.png',plotdir);
print(plotname,'-dpng','-r300');
fprintf('save plot to %s\n',plotname);

% KSN1 prl plot: contours in oscillation parameter space
% Lisa, October 2020

% work flow for SterileAnalysis clas:
% 1. set up a (Multi-)RunAnalysis object "R" with general settings (Runlist, FSD,...)
% 2. set up SterileAnalysis object "S" with "R" as input argument
% 3. from here you can load chi2 grids / plot contours, change some basic settings
%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType','Real',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

R = MultiRunAnalysis(RunAnaArg{:});
R.chi2 = 'chi2CMShape';
%% 2. SterileAnalysis class
SterileArg = {'RunAnaObj',R,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',40};

S = SterileAnalysis(SterileArg{:});

%% 3. contour plot in oscillation parameter space (you can switch on/off foreign contours, all on by default
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline';           % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some

Arg = {'SavePlot','OFF','BestFit','OFF','FinalSensitivity','ON','Style','PRL'};
%
mbbexplim = 0.165;
nGrid = 50;
nSamples = 5e3;
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_NuBetaBeta/results/'];
savenameNH = sprintf('%sksn1_NuBetaBeta_ToyMC_%.3feVlim_NH_Grid%.0f_Samples%.0f_New_light.mat',savedir,mbbexplim,nGrid,nSamples);
savenameIH = sprintf('%sksn1_NuBetaBeta_ToyMC_%.3feVlim_IH_Grid%.0f_Samples%.0f_New_light.mat',savedir,mbbexplim,nGrid,nSamples);
dNH = importdata(savenameNH);
dIH = importdata(savenameIH);


% plot
figure('Units','normalized','Position',[0.1,0.1,0.382,0.66]);%0.618]);
[~,sin2T4SqNH] = S.Convert2Osci('sinT4',dNH.sint2Sq);
[~,sin2T4SqIH] = S.Convert2Osci('sinT4',dIH.sint2Sq);
              

[lIH,aIH] = boundedline(sin2T4SqIH,dIH.m4Sq,dIH.m4Sq_err);
lIH.delete;
aIH.FaceColor = rgb('Silver');
aIH.FaceAlpha = 0.2;

[lNH,aNH] = boundedline(sin2T4SqNH,dNH.m4Sq,dNH.m4Sq_err);
lNH.delete;
aNH.FaceColor = rgb('Silver');
aNH.FaceAlpha = 0.6;

[legHandle,legStr] = S.ContourPlotOsci(Arg{:},'HoldOn','ON');


legHandle = {legHandle{:},aNH,aIH};
%legStr = {legStr{:},sprintf('KamLAND-Zen NH 90%% C.L.'),sprintf('KamLAND-Zen IH 90%% C.L.')};
legStr = {legStr{:},sprintf('0\\nu\\beta\\beta NH 90%% C.L.'),sprintf('0\\nu\\beta\\beta IH 90%% C.L.')};

leg = legend([legHandle{:}],legStr{:},'EdgeColor','none','Location','northoutside');
%%            
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir)
if contains(savenameNH,'New')
    plotname = sprintf('%sNuBetaBeta_PRL_%.3feV_mbb3nu_Cosmo.pdf',plotdir,mbbexplim);
else
    plotname = sprintf('%sNuBetaBeta_PRL_%.3feV_mbb3nu_NonDeg.pdf',plotdir,mbbexplim);
end
%print(plotname,'-dpng','-r350');
export_fig(plotname);
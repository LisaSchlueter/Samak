% ksn2 calculate chi2 grid search
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;
freePar = 'E0 Norm Bkg';
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5}};

S = SterileAnalysis(SterileArg{:});
%% %S.PlotTwinData('BestFit','ON','SavePlot','ON')
%% twins
S.RunAnaObj.DataType = 'Twin';
S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','OFF');
S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
S.ContourPlot; close;

% chi2_T = S.chi2;
% sin2T4_contourT = S.sin2T4_contour;
% mNu4Sq_contourT = S.mNu4Sq_contour;
% 
% chi2 = S.chi2;
% chi2ref = S.chi2_ref;
% sin2T4 = S.sin2T4;
% mNu4Sq = S.mNu4Sq;
% sin2T4_contour = S.sin2T4_contour;
% mNu4Sq_contour = S.mNu4Sq_contour;
% 
% Twin_mNu4Sq = 92.7;
% Twin_sin2T4 = 0.024;
% S.Twin_mNu4Sq = Twin_mNu4Sq;
% S.Twin_sin2T4 = Twin_sin2T4;
% S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','OFF');
% S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
% S.ContourPlot('BestFit','ON'); close;
% 
% chi2_1 = S.chi2;
% chi2ref_1 = S.chi2_ref;
% sin2T4_1 = S.sin2T4;
% mNu4Sq_1 = S.mNu4Sq;
% sin2T4_contour_1 = S.sin2T4_contour;
% mNu4Sq_contour_1 = S.mNu4Sq_contour;
% sin2T4_bf = S.sin2T4_bf;
% mNu4Sq_bf = S.mNu4Sq_bf;
%%
% data
S.RunAnaObj.DataType = 'Real';
S.LoadGridFile('CheckSmallerN','ON',S.LoadGridArg{:},'ExtmNu4Sq','ON');
S.Interp1Grid('RecomputeFlag','ON','MinM4Sq',1,'MaxM4Sq',40^2);
S.ContourPlot('BestFit','ON'); close;
chi2_D = S.chi2-min(min(S.chi2));
sin2T4_contourD = S.sin2T4_contour;
mNu4Sq_contourD = S.mNu4Sq_contour;
sin2T4_bfD = S.sin2T4_bf;
mNu4Sq_bfD = S.mNu4Sq_bf;
%%
chi2Diff = chi2_D-chi2_T;
chi2Diff(abs(chi2Diff)>10) = NaN;
GetFigure;
surf(S.sin2T4,S.mNu4Sq,chi2Diff,'EdgeColor','none');
hold on;
pD = plot3(sin2T4_contourD,mNu4Sq_contourD,99.*ones(numel(mNu4Sq_contourD),1),'k-','LineWidth',2);
if  contains(freePar,'mNu')
    pDbf = plot3(sin2T4_bfD,mNu4Sq_bfD,99,'kx','LineWidth',2);
end
pT = plot3(sin2T4_contourT,mNu4Sq_contourT,99.*ones(numel(mNu4Sq_contourT),1),':','LineWidth',2.5,'Color',rgb('White'));
set(gca,'XScale','log');set(gca,'YScale','log');
PrettyFigureFormat('FontSize',22);
view(2); grid off;

c = colorbar;
c.Label.String = sprintf('\\Delta\\chi^2 (Data - Twin)');
c.Label.FontSize = get(gca,'FontSize');
leg = legend([pT,pD],'Twin sensitivity 95% C.L.','Data exclusion 95% C.L.','Location','southwest'); PrettyLegendFormat(leg);
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
if  contains(freePar,'mNu')
    leg.Title.String =(sprintf('{\\itm}_\\nu free'));
    xlim([5e-03 0.5]);
    ylim([3 1600]);
else
    leg.Title.String =(sprintf('{\\itm}_\\nu fix'));
    xlim([2e-03 0.5]);
    ylim([1 1600]);
end
leg.Title.FontWeight = 'normal';
leg.Title.FontSize = get(gca,'FontSize');
%%
print(gcf,[S.DefPlotName,'DataTwinGrid.png'],'-dpng','-r350');

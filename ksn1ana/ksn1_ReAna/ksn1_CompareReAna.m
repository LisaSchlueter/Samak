% ksn2 calculate chi2 grid search
%% settings that might change
nGridSteps = 30;
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'E0 Norm Bkg';

%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.064;
end
Real = MultiRunAnalysis('RunList','KNM1',...
    'chi2',chi2,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'SysBudget',24,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'AngularTFFlag','OFF',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','OFF',...
    'BKG_PtSlope',0);
Real.exclDataStart = Real.GetexclDataStart(range);
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

S = SterileAnalysis(SterileArg{:});

%% KSN-1 with PRL config:
Real.ModelObj.BKG_PtSlope = 0;
Real.AngularTFFlag = 'OFF';
Real.FSDFlag = 'Sibille0p5eV';
Real.BKG_PtSlope = 0;
Real.ELossFlag = 'KatrinT2';
S.RunAnaObj.SysBudget = 24;
S.nGridSteps = 50;
% load file and interpolate
S.LoadGridFile('ExtmNu4Sq','OFF','mNu4SqTestGrid','OFF');
S.Interp1Grid;
p1 = S.ContourPlot('BestFit','ON');
% find best fit
mNu4Sq_bf_1 = S.mNu4Sq_bf;
sin2T4_bf_1 = S.sin2T4_bf;
chi2min_1 = S.chi2_bf;
chi2_Null_1 = S.chi2_Null;
[DeltaChi2_1, SignificanceBF_1] = S.CompareBestFitNull;

% KSN-1 re-ana
Real.AngularTFFlag = 'ON';
Real.FSDFlag = 'KNM2_0p1eV';
Real.ModelObj.BKG_PtSlope = -2.2*1e-06;
Real.ELossFlag = 'KatrinT2A20';
S.RunAnaObj.SysBudget = 200;
S.nGridSteps = 30;
% load file and interpolate
S.LoadGridFile('ExtmNu4Sq','OFF','mNu4SqTestGrid',2);
S.Interp1Grid;
p2 = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.','BestFit','ON');
% find best fit
mNu4Sq_bf_2 = S.mNu4Sq_bf;
sin2T4_bf_2 = S.sin2T4_bf;
chi2min_2 = S.chi2_bf;
chi2_Null_2 = S.chi2_Null;
[DeltaChi2_2, SignificanceBF_2] = S.CompareBestFitNull;

% legend
leg = legend([p1,p2],'PRL-2021','Re-Analysis2021');
leg.Title.String = sprintf('KSN-1 model and systematics');
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);


ylim([1 40^2]);
xlim([6e-03 0.5]);

%% save
savename = [extractBefore(S.DefPlotName,'FSD'),sprintf('%s_ReAna.png',chi2)];
print(savename,'-dpng','-r350');
fprintf('save plot to %s \n',savename);


%% best fit display
fprintf('--------------------------------------\n');
fprintf('--------------------------------------\n');
fprintf('Best fit:     PRL --- Re-Ana \n');
fprintf('--------------------------------------\n');
fprintf('m4Sq (eV^2): %.1f --- %.1f  \n',mNu4Sq_bf_1,mNu4Sq_bf_2);
fprintf('sint4Sq:    %.3f --- %.3f   \n',sin2T4_bf_1,sin2T4_bf_2);
fprintf('chi2min:     %.1f --- %.1f   \n',chi2min_1,chi2min_2);
fprintf('chi2 NH:     %.1f --- %.1f   \n',chi2_Null_1,chi2_Null_2);
fprintf('Delta chi2:   %.1f --- %.1f   \n',DeltaChi2_1,DeltaChi2_2);
fprintf('--------------------------------------\n');
fprintf('--------------------------------------\n');


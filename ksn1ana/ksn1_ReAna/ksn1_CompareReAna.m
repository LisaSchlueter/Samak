% ksn1 re-analysis:
% compare PRL config with re-analysis result 
% contour plot
%% settings that might change
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
S.LoadGridArg = {'ExtmNu4Sq','OFF'};

% load file and interpolate
S.LoadGridFile(S.LoadGridArg{:});

if contains(freePar,'mNu')
    S.InterpMode = 'Mix';
    S.Interp1Grid('Maxm4Sq',19^2);
else
    S.InterpMode = 'spline';
    S.Interp1Grid;
end

[p1,~,p1bf] = S.ContourPlot('BestFit','ON');

% find best fit
S.FindBestFit('Mode','Imp');
mNu4Sq_bf_1 = S.mNu4Sq_bf;
sin2T4_bf_1 = S.sin2T4_bf;
chi2min_1 = S.chi2_bf;
chi2_Null_1 = S.chi2_Null;
[DeltaChi2_1, SignificanceBF_1] = S.CompareBestFitNull;

%% KSN-1 re-ana
Real.AngularTFFlag = 'ON';
Real.FSDFlag = 'KNM2_0p1eV';
Real.ModelObj.BKG_PtSlope = -2.2*1e-06;
Real.ELossFlag = 'KatrinT2A20';
S.RunAnaObj.SysBudget = 200;
S.nGridSteps = 50;


% load file and interpolate
if contains(freePar,'mNu')
    S.LoadGridArg = {'mNu4SqTestGrid',2};
    S.InterpMode = 'Mix';
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',36^2);
else
    S.LoadGridArg = {'mNu4SqTestGrid',2};
    S.InterpMode = 'spline';
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid;
end


[p2,~,p2bf] = S.ContourPlot('HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.','BestFit','ON');

% find best fit
S.FindBestFit('Mode','Imp');
mNu4Sq_bf_2 = S.mNu4Sq_bf;
sin2T4_bf_2 = S.sin2T4_bf;
chi2min_2 = S.chi2_bf;
chi2_Null_2 = S.chi2_Null;
[DeltaChi2_2, SignificanceBF_2] = S.CompareBestFitNull;

ylim([1 40^2]);
xlim([6e-03 0.5]);

SigmaBF1 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_1,'nPar',2);
SigmaBF2 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_2,'nPar',2);
%% legend
ExtLeg = 'OFF';
if strcmp(ExtLeg,'ON')
    leg = legend([p1,p2,p1bf,p2bf],'PRL-2021','Re-Analysis2021',...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sq_bf_1,sin2T4_bf_1,SigmaBF1),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sq_bf_2,sin2T4_bf_2,SigmaBF2));
    leg.NumColumns = 2;
    extraStr = '_ExtLeg';
else
    leg = legend([p1,p2],'PRL-2021','Re-Analysis2021');
     extraStr = '';
end
leg.Title.String = sprintf('KSN-1 model and systematics');
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);


%% save
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_ReAna/plots/'];
MakeDir(savedir);
savename = [savedir,sprintf('ksn1_ReAna_ComparisonPlt_%s_%s_%.0feV.png',strrep(freePar,' ',''),chi2,range)];
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


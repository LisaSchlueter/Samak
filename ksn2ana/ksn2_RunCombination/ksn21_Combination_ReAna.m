% combine ksn1 and ksn2, nu-mass fixed
% simply add chi2-map2
% with KSN-1 Re-Ana chi2 map
%% settings that might change
chi2Name = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 50;
range = 40;
BF = 'ON';
% savedir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombinations/results'];
% savefile = sprintf('%sksn21_Combination_ReAna_%s_%s.mat',DataType,chi2);

%% configure RunAnalysis object
if strcmp(chi2Name,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2Name,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2Name,...
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
    'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};

%%
S = SterileAnalysis(SterileArg{:});
%% load ksn1
S.RunAnaObj.ModelObj.BKG_PtSlope = -2.2*1e-06;
S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 200;
S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
S.RunAnaObj.NonPoissonScaleFactor = 1.064;
S.RunAnaObj.chi2 = 'chi2CMShape';
S.nGridSteps = 50;
S.LoadGridArg = {'CheckLarger','ON','ExtmNu4Sq','OFF','mNu4SqTestGrid',2};

% load grid and plot
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('RecomputeFlag','ON');%,'Maxm4Sq',34.2^2);
[p1tot,~,p1bf] = S.ContourPlot('BestFit',BF,'CL',95,'HoldOn','OFF','Color',rgb('FireBrick'),'LineStyle',':');

% store best fit and contour
S.FindBestFit;
S.FindBestFit('Mode','Imp')
chi2_ref_1 = S.chi2_ref;
chi2_null_1 = S.chi2_Null;
mNuSq_bf_1 = S.mNuSq_bf;
sin2T4_bf_1 = S.sin2T4_bf;
mNu4Sq_bf_1 = S.mNu4Sq_bf;
% find significance
[~, SignificanceBF_1] = S.CompareBestFitNull;

% load again 
% restrict binning to match KSN2
S.LoadGridFile(S.LoadGridArg{:}); 
S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',38^2);
mNu4Sq_k1 = S.mNu4Sq;
sin2T4_k1 = S.sin2T4;
chi2_k1   = S.chi2;
chi2ref_k1= S.chi2_ref;
sum(sum(isnan(S.chi2)))
dof1 = S.dof;

%% load ksn2
S.nGridSteps = 30;
S.RunAnaObj.ModelObj.BKG_PtSlope = 3*1e-06;
S.RunAnaObj.FSDFlag = 'KNM2_0p5eV';
S.RunAnaObj.DataSet = 'Knm2';
S.RunAnaObj.RunData.RunName = 'KNM2_Prompt';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 40;
S.RunAnaObj.NonPoissonScaleFactor = 1.112;
S.RunAnaObj.chi2 = 'chi2CMShape';
LoadGridArg_i = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
S.LoadGridArg = {LoadGridArg_i{:},'Extsin2T4','ON'};

% load grid and plot
S.LoadGridFile(S.LoadGridArg{:});
S.InterpMode = 'Mix';
S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',36^2);
[p2tot,sinMin,p2bf] = S.ContourPlot('BestFit',BF,'CL',95,'HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.');

% store best fit and contour
S.FindBestFit;
S.FindBestFit('Mode','Imp')
chi2_ref_2 = S.chi2_ref;
chi2_null_2 = S.chi2_Null;
mNuSq_bf_2 = S.mNuSq_bf;
sin2T4_bf_2 = S.sin2T4_bf;
mNu4Sq_bf_2 = S.mNu4Sq_bf;
mNu4Sq_k2_contour = S.mNu4Sq_contour;
sin2T4_k2_contour = S.sin2T4_contour;
% find significance
[~, SignificanceBF_2] = S.CompareBestFitNull;

% load again with same binning as KSN-1 (restricted) to sum chi^2 maps
S.LoadGridArg = {LoadGridArg_i{:},'Extsin2T4','OFF'};
S.LoadGridFile(S.LoadGridArg{:}); 
S.InterpMode = 'spline';
S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',38^2);
mNu4Sq_k2 = S.mNu4Sq;
sin2T4_k2 = S.sin2T4;
chi2_k2   = S.chi2;
chi2ref_k2= chi2_ref_2;
dof2 = S.dof;
%% check if they have same binning
 if sum(sum(sin2T4_k2~=sin2T4_k1))>0 || sum(sum(mNu4Sq_k1~=mNu4Sq_k2))>0
     fprintf(2, 'KSN-1 and KSN-2 do not have  the same binning - return \n');
     return;
 end
 
%% combi 
% sum
chi2_sum = chi2_k1+chi2_k2;
S.chi2 = chi2_sum;
S.chi2_ref = min(min(chi2_sum));%chi2ref_k1+chi2ref_k2;
[p12tot,sinMin,p12bf] = S.ContourPlot('BestFit','ON','CL',95,'HoldOn','ON','Color',rgb('DodgerBlue'),'LineStyle','-');

% find best fit;
S.FindBestFit;
chi2_ref_12 = S.chi2_ref;
chi2_null_12 =  chi2_null_1+chi2_null_2;
mNuSq_bf_12 = S.mNuSq_bf;
sin2T4_bf_12 = S.sin2T4_bf;
mNu4Sq_bf_12 = S.mNu4Sq_bf;
S.chi2_Null = chi2_null_12;
% find significance
[~, SignificanceBF_12] = S.CompareBestFitNull;
 %% plot best fits on top (KSN-1 is hidden otherwise)
pbf1 = plot(sin2T4_bf_1,mNu4Sq_bf_1,'x','MarkerSize',9,'Color',p1tot.Color,'LineWidth',p1tot.LineWidth);

%% legend
SigmaBF1 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_1,'nPar',2);
SigmaBF2 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_2,'nPar',2);
SigmaBF12 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_12,'nPar',2);
%
ExtLeg = 'ON';
if strcmp(ExtLeg,'ON')
    leg = legend([p1tot,p2tot,p12tot,p1bf,p2bf,p12bf],'KSN-1','KSN-2','KSN-1 and KSN-2',...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sq_bf_1,sin2T4_bf_1,SigmaBF1),...
        sprintf('{\\itm}_4^2 = %.2f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sq_bf_2,sin2T4_bf_2,SigmaBF2),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sq_bf_12,sin2T4_bf_12,SigmaBF12));
    leg.NumColumns = 2;
    extraStr = '_ExtLeg';
else
    leg = legend([p1tot,p2tot,p12tot],'KSN-1','KSN-2','KSN-1 and KSN-2 combined');
    extraStr = '';
end
PrettyLegendFormat(leg);
legend boxoff

%% save
title(S.GetPlotTitle,'FontWeight','normal')
xlim([2e-03,0.5]);

 %% save plot
 savedir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/plots/'];
 MakeDir(savedir);
 savefile = sprintf('%sksn21_Combination_ReAna%s',savedir,extraStr);
 print([savefile,'.png'],'-dpng','-r350');
 
 ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
 export_fig([savefile,'.pdf']);
 fprintf('save plot to %s \n',savefile);
 %% save mini result file for compatibility test
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir)
savefile = sprintf('%sksn21_Combination_ReAna.mat',savedir);

dof12 = dof1+dof2+2;
save(savefile,'chi2_ref_1','chi2_ref_2','chi2_ref_12','dof1','dof2','dof12',...
    'chi2_null_1','chi2_null_2','chi2_null_12',...
    'mNu4Sq_bf_1','mNu4Sq_bf_2','mNu4Sq_bf_12',...
    'sin2T4_bf_1','sin2T4_bf_2','sin2T4_bf_12');



% combine ksn1 and ksn2, nu-mass free
% plot with ksn1, ksn2, ksn1+2
chi2          = 'chi2CMShape';
DataType      = 'Real';
nGridStepsCommon    = 30;
freePar       = 'mNu E0 Bkg Norm';
range         = 40;
savedir = [getenv('SamakPath'),'/SterileAnalysis/GridSearchFiles/Combi/',DataType,'/'];
MakeDir(savedir)
savename = sprintf('%sKSN12Combi_ReAna_GridSearch_%s_%s_Uniform_%s_%.0fnGrid.mat',...
    savedir,DataType,strrep(freePar,' ',''),chi2,nGridStepsCommon);

% load combi file
if exist(savename,'file') 
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
     fprintf('load results from file %s \n',savename)
    return
end

savename2 = sprintf('%sKSN12Combi_ReAna_RunAnaObj_%s.mat',savedir,DataType);
if exist(savename2,'file')
    A = importdata(savename2);
else
    % int obejct for interpolatio and plot
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    save(savename2,'A');
end

A.fixPar = ConvertFixPar('freePar','mNu E0 Norm Bkg','nPar',A.nPar,'nPixel',1);
A.chi2 = chi2;

%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
%% ksn1
KSN1config = 'New';
S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.NonPoissonScaleFactor = 1.064;
S.RunAnaObj.chi2 = 'chi2CMShape';
switch KSN1config
    case 'Old'
        % OLD settings (PRL publication)
        S.nGridSteps = 50;
        S.RunAnaObj.ELossFlag  = 'KatrinT2';
        S.RunAnaObj.AngularTFFlag ='OFF';
        S.RunAnaObj.SysBudget = 24;
        S.RunAnaObj.FSDFlag = 'Sibille0p5eV';
        S.RunAnaObj.ModelObj.BKG_PtSlope = 0;
        S.LoadGridArg = '';
        S.InterpMode = 'Mix';
        % load map
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',19^2);
    case 'New'
        % Re-Analysis
        S.nGridSteps = 50;
        S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
        S.RunAnaObj.AngularTFFlag ='ON';
        S.RunAnaObj.SysBudget = 200;
        S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
        S.RunAnaObj.ModelObj.BKG_PtSlope = -2.2*1e-06;
        S.LoadGridArg = {'mNu4SqTestGrid',2};
        S.InterpMode = 'Mix';
        % load map
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',36^2);
end


% store
mNu4Sq_k1 = S.mNu4Sq;
sin2T4_k1 = S.sin2T4;
chi2_k1   = S.chi2;
chi2ref_k1= S.chi2_ref;
sum(sum(isnan(S.chi2)))

% plot
[p1tot,~,pbf1] = S.ContourPlot('BestFit','ON','CL',95,'HoldOn','OFF','Color',rgb('FireBrick'),'LineStyle',':');
mNu4Sq_1_contour = S.mNu4Sq_contour;
sin2T4_k1_contour = S.sin2T4_contour;
dof1 = S.dof;

% find significance
[DeltaChi2_1, SignificanceBF_1] = S.CompareBestFitNull;
mNu4Sqbf_k1 = S.mNu4Sq_bf;
sin2T4bf_k1 = S.sin2T4_bf;
%% load ksn2
S.InterpMode = 'spline';
S.nGridSteps = 50;
S.RunAnaObj.FSDFlag = 'KNM2_0p5eV';
S.RunAnaObj.DataSet = 'Knm2';
S.RunAnaObj.RunData.RunName = 'KNM2_Prompt';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 40;
S.RunAnaObj.NonPoissonScaleFactor = 1.112;
S.RunAnaObj.ModelObj.BKG_PtSlope = 3.*1e-06;
% stat and syst
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridFile('CheckLarger','OFF','mNu4SqTestGrid',5,'ExtmNu4Sq','OFF');            
S.Interp1Grid('RecomputeFlag','ON');
mNu4Sq_k2 = S.mNu4Sq;
sin2T4_k2 = S.sin2T4;
chi2_k2   = S.chi2;
chi2ref_k2= S.chi2_ref;
mNu4Sq_k2_contour = S.mNu4Sq_contour;
sin2T4_k2_contour = S.sin2T4_contour;
%% get best fit from fine grid
savedir = sprintf('%sksn2ana/ksn2_BestFit/results/',getenv('SamakPath'));
savefile2 = sprintf('%sksn2_ImpBestFit_%s_%s_%s_nGridStepsFull%.0f_nGridStepsImp%.0f.mat',...
    savedir,DataType,'mNuE0NormBkg',chi2,40,25);
if exist(savefile2,'file') 
    d2 = importdata(savefile2);
    S.mNu4Sq_bf = d2.mNu4Sq_bf_inter;
    S.sin2T4_bf = d2.sin2T4_bf_inter;
    S.chi2_bf = d2.chi2_bf_inter;
else
    S.FindBestFit;
end
mNu4Sqbf_k2 = S.mNu4Sq_bf;
sin2T4bf_k2 = S.sin2T4_bf;

[p2tot,sinMin,pbf2] = S.ContourPlot('BestFit','ON','CL',95,'HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.','ReCalcBF','OFF');
dof2 = S.dof;
% find significance
[DeltaChi2_2, SignificanceBF_2] = S.CompareBestFitNull;
%% KSN1+2
S.sin2T4 = d.sin2T4;
S.mNu4Sq = d.mnu4Sq;
S.chi2 = d.chi2;
S.chi2_ref = d.chi2_ref;
S.chi2_Null = d.FitResult_Null.chi2min;
S.mNuSq = cell2mat(cellfun(@(x) x.par(1),d.FitResults,'UniformOutput',0));
S.E0 = cell2mat(cellfun(@(x) x.par(2),d.FitResults,'UniformOutput',0));
dof12 = 46;
S.dof = dof12;
S.InterpMode = 'spline';
S.Interp1Grid('Maxm4Sq',40^2);%,'Minm4Sq',40);
[p12tot,~,pbf12] = S.ContourPlot('HoldOn','ON','BestFit','ON','SavePlot','OFF');

mNu4Sq_k12_contour = S.mNu4Sq_contour;
sin2T4_k12_contour = S.sin2T4_contour;

mNu4Sqbf_k12 = S.mNu4Sq_bf;
sin2T4bf_k12 = S.sin2T4_bf;

% find significance
[DeltaChi2_tot, SignificanceBF_tot] = S.CompareBestFitNull;

% convert cl to sigma
SigmaBF1 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_1,'nPar',2);
SigmaBF2 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_2,'nPar',2);
SigmaBF12 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_tot,'nPar',2);
%%
ExtLeg = 'ON';
legStr = sprintf('%s: {\\itm}_\\nu^2 free at %.0f%% C.L.',S.GetPlotTitle('Mode','data'),S.ConfLevel);
if strcmp(ExtLeg,'ON')
    leg = legend([p1tot,p2tot,p12tot,pbf1,pbf2,pbf12],...
        'KSN-1','KSN-2','KSN-1 and KSN-2',...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k1,sin2T4bf_k1,SigmaBF1),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k2,sin2T4bf_k2,SigmaBF2),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k12,sin2T4bf_k12,SigmaBF12));
    leg.NumColumns = 2;
    title(legStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
else
    leg = legend([p1tot,p2tot,p12tot],'KSN-1','KSN-2','KSN-1 and KSN-2 combined');
    leg.Title.String = legStr;
    leg.Title.FontWeight = 'normal';
    title('');
end
PrettyLegendFormat(leg);

xlim([2e-03,0.5]);
if strcmp(ExtLeg,'ON')
    ylim([0.5 42^2]);
else
    ylim([1 40^2]);
end
%%
% save plot
plotname = [extractBefore(S.DefPlotName,'Grid'),sprintf('KSN12_%s_Combination_mNuE0NormBkg_%s.png',...
    S.GetPlotTitle('Mode','data'),chi2)];
if strcmp(ExtLeg,'ON')
    plotname = strrep(plotname,'.png','_BF.png');
end
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname)
 ylabel(sprintf('{\\itm}_4^2 (eV^{ 2})'));
export_fig(gcf,strrep(plotname,'.png','.pdf'));

%% save combi file
savedir = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
MakeDir(savedir);
savefile = sprintf('%sksn21_Combi_freemNuSq_ReAna.mat',savedir);
save(savefile,'chi2ref_k1','chi2ref_k2','DeltaChi2_1','DeltaChi2_2','DeltaChi2_tot',...
    'mNu4Sq_1_contour','mNu4Sq_k2_contour','sin2T4_k1_contour','sin2T4_k2_contour',...
    'mNu4Sq_k12_contour','sin2T4_k12_contour','d',...
    'mNu4Sqbf_k1','mNu4Sqbf_k2','mNu4Sqbf_k12',...
    'sin2T4bf_k1','sin2T4bf_k2','sin2T4bf_k12',...
    'dof1','dof2','dof12');

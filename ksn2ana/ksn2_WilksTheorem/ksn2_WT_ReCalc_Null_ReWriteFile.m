% recalculate chi2 for MC truth were re-calculated
% look and results and overwrite file

Hypothesis = 'H1';
randMC = [1:139,577:757];%130;%[1:201,534:748]; % still need to be re-calc: 130-139,534:748
Twin_sin2T4 = 0.0240;
Twin_mNu4Sq = 92.7;
chi2 = 'chi2CMShape';
DataType    = 'Twin';
freePar     = 'E0 Norm Bkg';
nGridSteps  = 25;
range       = 40;

%% get randomized MC data
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',...
    savedir,chi2,range,Twin_mNu4Sq,Twin_sin2T4,5000);

if exist(savefile,'file')
    d = importdata(savefile);
else
    fprintf('file %s not found \n',savefile)
    return
end

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
A.exclDataStart = A.GetexclDataStart(range);
TBDIS_i = A.RunData.TBDIS;
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC','OFF',...
    'Twin_mNu4Sq',Twin_mNu4Sq,...
    'Twin_sin2T4',Twin_sin2T4,...
    'InterpMode','lin'};

S = SterileAnalysis(SterileArg{:});
%%

ReCalc_chi2Null = zeros(numel(randMC),1);
chi2Null        = zeros(numel(randMC),1);

progressbar('rewrite files')
for i=1:numel(randMC)
    progressbar(i/numel(randMC));
    S.RandMC = randMC(i);
    S.RandMC_TBDIS =   d.TBDIS_mc(:,randMC(i));
    filename = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');

    dtmp = importdata(filename);
    ReCalc_chi2Null(i) = dtmp.ReCalc_chi2Null;
    chi2Null(i) = dtmp.FitResults_Null.chi2min;
    
    % change file
    FitResults_Null = dtmp.FitResults_Null;
    FitResults_Null.chi2min = ReCalc_chi2Null(i);
    chi2_H0 = chi2Null(i);
    
    save(filename,'FitResults_Null','chi2_H0','-append');
end


%%
dof = dtmp.FitResults{1}.dof;
GetFigure;
histogram(chi2Null,'Normalization','probability','BinWidth',3);
hold on;
hchi2 = histogram(ReCalc_chi2Null,'Normalization','probability','BinWidth',3);
x = linspace(0,dof*3,1e3);
y = chi2pdf(x,dof);
pchi2 = plot(x,y*hchi2.BinWidth,'Color',rgb('Black'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\chi^2_{null} (%.0f dof)',dof));
ylabel('Frequency');
 
 


% ksn2 proposed model blinding
% ksn2 calculate chi2 grid search with different settings
%% settings that might change
Mode                  = 'Compute'; %'Compute';
chi2                  = 'chi2CMShape';
DataType              = 'Twin';
nGridSteps            = 25;
range                 = 40;
NonPoissonScaleFactor = 1.112;
Twin_mNu4Sq           = 10^2;
Twin_sin2T4           = 0.03;
SysBudget             = [40,446];%,441,442,443];
RandMC                = 1;
%% configure RunAnalysis object
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
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

%% calculate random twin used for all budgets
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_ModelBlinding/results/'];
savename = sprintf('%sksn2_RandMCTBDIS_mNu4Sq%.3g_sinT4Sq%.3g_%s_SysBudget%.0f.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,chi2,A.SysBudget);
if exist(savename,'file')
    d = importdata(savename);
    TBDIS_mc = d.TBDIS_mc;
else
    A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    
    A.ModelObj.BKG_RateSec_i = A.ModelObj.BKG_RateSec;
    A.ModelObj.normFit_i = A.ModelObj.normFit;
    if Twin_mNu4Sq~=0 || Twin_sin2T4~=0
        A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);
    end
    A.ModelObj.ComputeTBDDS;
    A.ModelObj.ComputeTBDIS;
    TBDIS_i = A.ModelObj.TBDIS';
    
    TBDIS_mc = mvnrnd(TBDIS_i,A.FitCMShape,1)';
    FitCMShape = A.FitCMShape;
    qU = A.ModelObj.qU;
    save(savename,'TBDIS_mc','FitCMShape','qU')
end
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC',RandMC,...
    'RandMC_TBDIS',TBDIS_mc,...
    'Twin_mNu4Sq',Twin_mNu4Sq,...
    'Twin_sin2T4',Twin_sin2T4,...
    };

%%
S = SterileAnalysis(SterileArg{:});
HoldOn='OFF';
pHandle = cell(numel(SysBudget),1);
legStr  = cell(numel(SysBudget),1);
Colors = jet(numel(SysBudget));
mNu4Sq_bf = zeros(numel(SysBudget),1);
sin2T4_bf = zeros(numel(SysBudget),1);

LineStyles = {'-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
for i=1:numel(SysBudget)
    A.SysBudget = SysBudget(i);
    switch Mode
        case 'Compute'
            A.ComputeCM;
            S.GridSearch;
        case 'Display'
            S.LoadGridFile;
            S.Interp1Grid('Maxm4Sq',34.^2,'nInter',1e3);
            pHandle{i} = S.ContourPlot('HoldOn',HoldOn,'BestFit','ON','CL',95,'Color',Colors(i,:),'LineStyle',LineStyles{i});
            HoldOn='ON';
             S.FindBestFit('Mode','Def')
             S.FindBestFit('Mode','Imp')
            sin2T4_bf(i) = S.sin2T4_bf;
            mNu4Sq_bf(i) = S.mNu4Sq_bf;
            if SysBudget(i)==40
                legStr{i} = sprintf('No additional uncertainty');
            elseif  SysBudget(i)==440
                legStr{i} = sprintf('\\mu = 0.05 %% \\sigma = 0.1 %%');
            elseif  SysBudget(i)==441
                legStr{i} = sprintf('\\mu = 0.1 %% \\sigma = 0.1 %%');
            elseif  SysBudget(i)==442
                legStr{i} = sprintf('\\mu = 0.2 %% \\sigma = 0.1 %%');
            elseif  SysBudget(i)==443
                legStr{i} = sprintf('\\mu = 0.4 %% \\sigma = 0.1 %%');
            end
    end
end

if strcmp(Mode,'Display')
    leg = legend([pHandle{:}],legStr);
    PrettyLegendFormat(leg);
    leg.Title.String = sprintf('Add. random uncorrelated uncertainty');
    leg.Title.FontWeight = 'normal';
    title(sprintf('MC truth: {\\itm}_4^2 = %.3geV^2 , |{\\itU}_{e4}|^2 = %.3g',Twin_mNu4Sq,Twin_sin2T4),...
        'FontWeight','normal','FontSize',get(gca,'FontSize'));
end

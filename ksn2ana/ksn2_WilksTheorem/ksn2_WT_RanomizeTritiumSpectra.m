%  generate random tritium spectra locally and give as input to grid search
%  in Wilk's theorem test
NrandMC = 5000;
chi2 = 'chi2CMShape';
DataType = 'Twin';
range = 40;

Hypothesis = 'H1';
switch Hypothesis
    case 'H0'
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
    case 'H1'
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_NullHypothesis_%.0fsamples.mat',savedir,chi2,range,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',...
        savedir,chi2,range,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file')
    load(savefile);
else
    %% configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
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
    
    % set range
    A.exclDataStart = A.GetexclDataStart(range);
    
    % get covariance matrix
    FitResults_i =  A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    nPixels = A.ModelObj.nPixels;
    Par_i = FitResults_i.par;
    
    
    if strcmp(A.chi2,'chi2CMShape')
        A.ComputeCM('BkgCM','ON','BkgPtCM','ON');
    end
    
    %% get expected TBDIS_i
    if Twin_mNu4Sq~=0 || Twin_sin2T4~=0
        A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);
        A.ModelObj.ComputeTBDDS(...
            'E0_bias',Par_i(2),...
            'B_bias',Par_i(3:2+nPixels),...
            'N_bias',Par_i(3+nPixels:3+2*nPixels-1));
        A.ModelObj.ComputeTBDIS;
        TBDIS_i = A.ModelObj.TBDIS(A.exclDataStart:end)';
    else
        TBDIS_i = A.ModelObj.TBDIS(A.exclDataStart:end)';
    end
    
    %% randomize spectra
    FitCMShape =  A.FitCMShape(A.exclDataStart:end,A.exclDataStart:end);
    qU = A.ModelObj.qU(A.exclDataStart:end);
    
    TBDIS_mc = zeros(A.ModelObj.nqU,NrandMC);
    TBDIS_mc(A.exclDataStart:end,:) =  mvnrnd(TBDIS_i,FitCMShape,NrandMC)';
    
    exclDataStart = A.exclDataStart;
    %% save
    save(savefile,'TBDIS_mc','TBDIS_i','FitCMShape','qU','RunAnaArg','exclDataStart')
end

%% sanity check
GetFigure;
[a,b] = boundedline(qU-18574,mean(TBDIS_mc(exclDataStart:end,:)')',50.*std(TBDIS_mc(exclDataStart:end,:)'));
hold on
e1 = errorbar(qU-18574,TBDIS_i,50.*sqrt(TBDIS_i), '.','CapSize',0,'LineWidth',1.5);

xlim([-44 50])
PrettyFigureFormat;
xlabel('Retarding potential -18574 (eV)');
ylabel('Counts');

leg = legend([e1,b],'Asimov model',sprintf('1\\sigma error band'));
PrettyLegendFormat(leg);

%% check if there are duplicates
   [TBDIS_mcNonId,IdxNonId] = unique(TBDIS_mc','rows');



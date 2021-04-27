% look at randomized spectra and find out whats wrong
Hypothesis = 'H0';
switch Hypothesis
    case 'H0'
        randMC = 1:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
        randMC = [1:340,384:794];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end
DataType = 'Twin';
freePar = 'E0 Norm Bkg';
nGridSteps = 25;
range = 40;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_RandSpectra_NullHypothesis_%.0fsamples.mat',savedir,numel(randMC));
else
    savefile = sprintf('%sksn2_WilksTheorem_RandSpectra_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,numel(randMC));
end

if exist(savefile,'file')
    d = importdata(savefile);
else
    savedirGrids = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm2/TwinRandomizedMC/'];
    
    savenameCommon =   [savedirGrids,'KSN2_GridSearch_KNM2_Prompt_Twin_E0BkgNorm_40eVrange_FSDKNM2_0p5eV_',...
        chi2,'_25nGrid_KatrinT2A20_AngTF'];
    if strcmp(Hypothesis,'H1')
        savenameCommon = [savenameCommon,sprintf('_sinT4Sq%.3f_mNu4Sq%.1f',Twin_sin2T4,Twin_mNu4Sq)];
    end
    if strcmp(chi2,'chi2CMShape')
        savenameCommon = strrep(savenameCommon,'25nGrid','25nGrid_Budget40');
    end
    
    TBDIS_mc = zeros(numel(randMC),38);
    for i=1:numel(randMC)
        progressbar(i/numel(randMC));
        df = [savenameCommon,'_RandMC',num2str(randMC(i)),'_ExtmNu4Sq_mNu4SqTestGrid2.mat'];
        d = importdata(df);
        TBDIS_mc(i,:) = d.TBDIS_mc;
    end
    save(savefile,'TBDIS_mc')
end

if ~isfield(d,'TBDIS_Asimov')
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
    A.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    if Twin_mNu4Sq~=0 || Twin_sin2T4~=0
        A.ModelObj.BKG_RateSec_i = A.ModelObj.BKG_RateSec;
        A.ModelObj.normFit_i = A.ModelObj.normFit;
        A.ModelObj.SetFitBiasSterile(Twin_mNu4Sq,Twin_sin2T4);
        A.ModelObj.ComputeTBDDS;
        A.ModelObj.ComputeTBDIS;
        TBDIS_Asimov = A.ModelObj.TBDIS';
    else
        %                         A.ModelObj.ComputeTBDDS;
        %                         A.ModelObj.ComputeTBDIS;
        TBDIS_Asimov = A.ModelObj.TBDIS';
    end
    
    qU = A.ModelObj.qU;
    Time = A.ModelObj.qUfrac.*A.ModelObj.TimeSec;
    save(savefile,'TBDIS_Asimov','qU','Time','-append')
    d = importdata(savefile);
end
%% display
TBDIS_Asimov = d.TBDIS_Asimov;
TBDIS_mc = d.TBDIS_mc;
exclIdx = 11;

GetFigure;
[a,b] = boundedline(d.qU(exclIdx:end)-18574,mean(d.TBDIS_mc(:,exclIdx:end)),50.*std(d.TBDIS_mc(:,exclIdx:end)));
hold on
e1 = errorbar(d.qU(exclIdx:end)-18574,d.TBDIS_Asimov(exclIdx:end),50.*sqrt(d.TBDIS_Asimov(exclIdx:end)), '.','CapSize',0,'LineWidth',1.5);

xlim([-44 50])
PrettyFigureFormat;
xlabel('Retarding potential -18574 (eV)');
ylabel('Counts');
%%
SumCounts    = sum(d.TBDIS_Asimov(exclIdx:end));
SumCounts_mc = sum(d.TBDIS_mc(:,exclIdx:end),2);

GetFigure;
histogram(SumCounts_mc)

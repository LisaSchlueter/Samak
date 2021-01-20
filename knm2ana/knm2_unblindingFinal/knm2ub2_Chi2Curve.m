range        = 40;
PlotChi2Scan = 'ON';
chi2         = 'chi2Stat';%CMShape';
freePar      = 'mNu E0 Bkg Norm qU';
SysBudget     = 38;
DataType      = 'Real';
AnaFlag       = 'Ring';%StackPixel';

if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/Chi2Curve/'];

savename = sprintf('%sknm2ub2_Chi2Curve_%s_%.0feV_%s_%s_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2,AnaFlag);
if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if exist(savename,'file') && 1==2
    load(savename,'ScanResult','FitResult','A')
    fprintf('load %s  \n',savename);
else
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'RingMerge','Full',...
        'pullFlag',99,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq)};

    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    %% Chi2 - scan
    A.Fit;
    FitResult = A.FitResult;
    ScanResult = A.GetAsymFitError('Mode','Uniform');
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A','SigmaSq','ScanResult')
end

if strcmp(PlotChi2Scan,'ON')
    A.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult,'FitResult',A.FitResult);
end

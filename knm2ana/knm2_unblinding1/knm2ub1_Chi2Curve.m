range        = 40;
PlotChi2Scan = 'ON';
chi2         = 'chi2CMShape';%Stat';
freePar      = 'mNu E0 Bkg Norm';
SysBudget     = 38;
DataType      = 'Twin';
AnaFlag       = 'StackPixel';

savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/'];

savename = sprintf('%sknm2ub1_Chi2Curve_%s_%.0feV_%s_%s_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),chi2,AnaFlag);
if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%.0f.mat',SysBudget));
end

if exist(savename,'file') 
    load(savename,'ScanResult','FitResult','A')
else
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'RingMerge','Full',...
        'chi2',chi2,...
        'pullFlag',99,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq)};
    
    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    if ~strcmp(chi2,'chi2Stat')
        A.NonPoissonScaleFactor = 1.112;
        A.SetNPfactor; % convert to right dimension (if multiring)
        A.chi2 = chi2;
        A.ComputeCM;
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
    % A.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult2,'FitResult',A.FitResult);
end

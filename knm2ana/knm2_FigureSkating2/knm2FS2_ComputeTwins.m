% KNM2 Figure skating twins
range = 40;
Chi2Scan = 'ON';

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');

RunAnaArg = {'RunList','KNM2_Prompt',...     % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...             % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...            % energy-loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...               % FPD segmentations -> pixel combination
    'chi2','chi2CMShape',...                    % statistics only
    'NonPoissonScaleFactor',1,...
    'TwinBias_Q',E0,...
    'ROIFlag','Default',...
    'fitter','minuit',...
    'SysBudget',36};

%% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
A.Fit;
A.PlotFit

%% Chi2 - scan
if strcmp(Chi2Scan,'ON')
    savedir = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating2/results/'];
    MakeDir(savedir);
    savename = sprintf('%sknm2FS2_ComputeTwinScanResults_%s_%.0feV.mat',savedir,A.chi2,range);
    
    if exist(savename,'file')
        load(savename,'ScanResult','ScanResult2','FitResult')
    else
        FitResult = A.FitResult;
        
        ScanResult = A.GetAsymFitError('Mode','Uniform');
        ScanResult2 = A.GetAsymFitError('Mode','Uniform','ParScanMax',0.1,'nFitMax',10);
        
        save(savename,'ScanResult','ScanResult2','FitResult')
    end
    
    A.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult,'FitResult',A.FitResult);
   % A.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult2,'FitResult',A.FitResult);
end

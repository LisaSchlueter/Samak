% KNM2 Figure skating twins

%MR.ComputeCM('SysEffect',struct('LongPlasma','ON'),'InitNormFit','OFF')
range = 40;
Chi2Scan = 'ON';

E0 = knm2FS_GetE0Twins('SanityPlot','OFF','Mode','FS2');
chi2 = 'chi2CMShape';

RunAnaArg = {'RunList','KNM2_Prompt',...
    'chi2','chi2Stat',...
    'DataType','Twin',...
    'fixPar','mNu E0 Bkg Norm',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'SysBudget',36,...
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2','chi2Stat',...
    'pullFlag',4,...
    'TwinBias_Q',E0,...
    'ROIFlag','Default',...
    'MosCorrFlag','OFF',...
    'NonPoissonScaleFactor',1};

%% build object of MultiRunAnalysis class
MR = MultiRunAnalysis(RunAnaArg{:});
MR.exclDataStart = MR.GetexclDataStart(range);

if ~strcmp(chi2,'chi2Stat')
    MR.NonPoissonScaleFactor = 1.112;
    MR.SetNPfactor; % convert to right dimension (if multiring)
    MR.chi2 = chi2;
    MR.ComputeCM;
end

%%
MR.Fit;
%MR.PlotFit

%% Chi2 - scan
if strcmp(Chi2Scan,'ON')
    savedir = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating2/results/'];
    MakeDir(savedir);
    savename = sprintf('%sknm2FS2_ComputeTwinScanResults_%s_%.0feV_MultiRing.mat',savedir,MR.chi2,range);
    
    if exist(savename,'file')
        load(savename,'ScanResult','FitResult')
    else
        FitResult = MR.FitResult;
        
        ScanResult = MR.GetAsymFitError('Mode','Uniform');
        %  ScanResult2 = MR.GetAsymFitError('Mode','Uniform','ParScanMax',0.1,'nFitMax',10);
        
        save(savename,'ScanResult','FitResult')
    end
    
    MR.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult,'FitResult',MR.FitResult);
    % A.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResult2,'FitResult',A.FitResult);
end

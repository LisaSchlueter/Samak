% KNM2 Figure skating twins

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    }; 

%% build object of MultiRunAnalysis class
D14 = MultiRunAnalysis(RunAnaArg{:},'ROIFlag','14keV');
D22 = MultiRunAnalysis(RunAnaArg{:},'ROIFlag','Default');

D14.exclDataStart = D14.GetexclDataStart(40);
D22.exclDataStart = D22.GetexclDataStart(40);
%%
FitResults14= D14.PlotFitRunList('Parameter','B','YLim',[160 280]);
FitResults22= D22.PlotFitRunList('Parameter','B','YLim',[160 280]);

%%
fprintf('Mean Background 14keV ROI = %.1f mcps \n',1e3.*mean(FitResults14.B));
fprintf('Mean Background 22keV ROI = %.1f mcps \n',1e3.*mean(FitResults22.B));

stairs(FitResults14.B-FitResults22.B);
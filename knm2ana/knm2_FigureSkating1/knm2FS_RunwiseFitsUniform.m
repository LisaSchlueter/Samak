% KNM2 Figure skating twins

range = 40;

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'ROIFlag','14keV'}; 

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});
D.exclDataStart = D.GetexclDataStart(range);

%%
FitResults= D.PlotFitRunList('Parameter','B','YLim',[170 280],'saveplot','pdf');
 D.PlotFitRunList('Parameter','E0','YLim',[-0.7,1],'DisplayStyle','Rel','saveplot','pdf');
D.PlotFitRunList('Parameter','N','YLim',[0.9,1.1],'saveplot','pdf');
D.PlotFitRunList('Parameter','pVal','YLim',[-0.2,1.2],'saveplot','pdf');
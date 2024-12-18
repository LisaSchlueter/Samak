% generate KNM2 twins with a plasma slope/steps

%% set up model to get start time stamps
[E0,~,~] = ConvertPlasmaDriftDay2E0Run('DriftPerDay',[6*1e-03, 0, 6*1e-03]',...
    'E0OffseteV',[0,0.1,-0.1]',...
    'E0ref',18573.70,...
    'SanityPlot','OFF');

RunArg = {'RunList','KNM2_Prompt',...
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'TwinBias_Q',E0,...
    'TwinBias_mNuSq',0};

MR = MultiRunAnalysis(RunArg{:});


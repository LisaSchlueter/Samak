% apply FSD correction to multiring fit
% idea: 3 periods with shifted rear wall/plasma potential + each period has a linear drift 

range = 40;
% read data and set up model
RunArg = {'RunList','KNM2_RW3',...
    'chi2','chi2Stat','DataType','Twin',...
    'fixPar','E0 Norm Bkg qU',...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2','chi2Stat',...
    'pullFlag',4,...
    'TwinBias_Q',18573.70};

R = MultiRunAnalysis(RunArg{:});
R.exclDataStart = R.GetexclDataStart(range);

%% broaden FSDs
MultiWeights = repmat(knm2_RWcombi_GetMultiWeights,[1,R.nRings]); % nPeaks x nRings
MultiPos     = repmat(knm2_RWcombi_GetMultiPos_E0fit,[1,R.nRings]); % nPeaks x nRings
TimeSec = MultiWeights.*R.RunData.TimeSec;
Slope = [10,20,30;12,22,32;14,24,34;16,26,36]'.*1e-03; % slope in eV/day
RectWidth = ConvertPlasmaDrift2Rect(Slope,TimeSec); % nPeaks x nRings
%%
FSDArg = {'MultiPos',MultiPos,'MultiWeights',MultiWeights,...
    'Sigma',RectWidth,'Dist','Rect','SanityPlot','OFF','ZoomPlot','OFF'};
R.ModelObj.LoadFSD(FSDArg{:});
%%
R.ModelObj.ComputeTBDDS;
R.ModelObj.ComputeTBDIS;
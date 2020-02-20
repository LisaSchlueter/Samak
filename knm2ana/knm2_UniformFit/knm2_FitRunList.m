function knm2_FitRunList(varargin)
% Fit all runs of KNM1
% Uniform 
% Feb 2020, Lisa

p = inputParser;
p.addParameter('RunList','KNM2_RW1',@(x)ischar(x));
p.addParameter('freePar','E0 Bkg Norm',@(x)ischar(x));
p.addParameter('Range',40,@(x)isfloat(x));              % fit range in eV below E0
p.addParameter('ROIFlag','Default',@(x)ismember(x,{'Default','14keV'}));   
p.parse(varargin{:});
RunList = p.Results.RunList;
freePar = p.Results.freePar;
Range   = p.Results.Range;
ROIFlag = p.Results.ROIFlag;

RunAnaArg = {'RunList',RunList,...  % define run number -> see GetRunList
    'fixPar',freePar,...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'ROIFlag',ROIFlag};

%% build object of MultiRunAnalysis class
A = MultiRunAnalysis(RunAnaArg{:});

%% modify some parameters in your analysis       
A.exclDataStart = A.GetexclDataStart(Range); % find correct data, where to cut spectrum

%% Fit
A.FitRunList;
end

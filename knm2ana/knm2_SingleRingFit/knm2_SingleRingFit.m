function R = knm2_SingleRingFit(varargin)
% single ring fit to KNM2 data
% analyze every rear wall period separately 
% 4 pseudo rings
p = inputParser;
p.addParameter('RunList','KNM2_RW1',@(x)ischar(x));
p.addParameter('freePar','mNu E0 Bkg Norm',@(x)ischar(x));
p.addParameter('Range',90,@(x)isfloat(x));              % fit range in eV below E0
p.addParameter('ROIFlag','14keV',@(x)ismember(x,{'Default','14keV'}));  
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat','chi2CMShape'}));  
p.addParameter('RingMerge','Full',@(x)ismember(x,{'Default','None','Full','Half','Azi'}));
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('MosCorrFlag','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

RunList   = p.Results.RunList;
freePar   = p.Results.freePar;
Range     = p.Results.Range;
ROIFlag      = p.Results.ROIFlag;
chi2         = p.Results.chi2;
RingMerge     = p.Results.RingMerge;
RecomputeFlag = p.Results.RecomputeFlag;
MosCorrFlag   = p.Results.MosCorrFlag;


%% settings
RunAnaArg = {'RunList',RunList,...  % define run number -> see GetRunList
    'fixPar',freePar,...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'RingMerge',RingMerge,...             % 'Full' == 4 Pseudo rings
    'fitter','minuit',...
    'minuitOpt','min;minos',...
    'ROIFlag',ROIFlag,...
    'chi2',chi2,...
    'MosCorrFlag',MosCorrFlag};             
            


%% read data and set up model: MultiRunAnalysis
A = MultiRunAnalysis(RunAnaArg{:}); % object of class MultiRunAnalysis
A.exclDataStart = A.GetexclDataStart(Range);

%% start ringwise analysis
R = RingAnalysis('RunAnaObj',A,'RingList',1:4); % object of class RingAnalysis

%% fit every ring - one after the other
R.FitRings('SaveResult','ON',...  
          'RecomputeFlag',RecomputeFlag,...  % load from storage or recalculate
          'AsymErr','OFF');         % asymmetric from scan and more correct uncertainties -> only for mNuSq

%% display
R.PlotFits('SavePlot','ON',...
          'Blind','ON',...       % show relative or absolute values
          'PlotPar',2,...        % 1 == neutrino mass, 2 == E0
          'YLim',[-0.2,0.15],... % force y-axis to good limits
          'linFit','ON');        % show linear fit
end

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

RunList       = p.Results.RunList;
freePar       = p.Results.freePar;
Range         = p.Results.Range;
ROIFlag       = p.Results.ROIFlag;
chi2          = p.Results.chi2;
RingMerge     = p.Results.RingMerge;
RecomputeFlag = p.Results.RecomputeFlag;
MosCorrFlag   = p.Results.MosCorrFlag;


%% settings
 SigmaSq =  0.0124+0.0025;
 if strcmp(chi2,'chi2Stat')
     NP = 1;
 else
     NP = 1.112;
 end
RunAnaArg = {'RunList',RunList,...  % define run number -> see GetRunList
    'fixPar',freePar,...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'RingMerge',RingMerge,...             % 'Full' == 4 Pseudo rings
    'fitter','minuit',...
    'minuitOpt','min;minos',...
    'chi2',chi2,...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','ON',...
    'SysBudget',38,...
    'FSD_Sigma',sqrt(SigmaSq),...
    'DopplerEffectFlag','FSD',...
    'NonPoissonScaleFactor',NP};             
            


%% read data and set up model: MultiRunAnalysis
A = MultiRunAnalysis(RunAnaArg{:}); % object of class MultiRunAnalysis
A.exclDataStart = A.GetexclDataStart(Range);

%% start ringwise analysis

R = RingAnalysis('RunAnaObj',A,'RingList',A.RingList); % object of class RingAnalysis

%% fit every ring - one after the other
R.FitRings('SaveResult','ON',...  
          'RecomputeFlag',RecomputeFlag,...  % load from storage or recalculate
          'AsymErr','OFF');                 % asymmetric from scan and more correct uncertainties -> only for mNuSq

%% display
  R.PlotFits('SavePlot','ON',...
          'Blind','ON',...       % show relative or absolute values
          'PlotPar',1,...        % 1 == neutrino mass, 2 == E0
          'YLim',[-3,3],... % force y-axis to good limits
          'linFit','ON');
      
R.PlotFits('SavePlot','ON',...
          'Blind','ON',...       % show relative or absolute values
          'PlotPar',2,...        % 1 == neutrino mass, 2 == E0
          'YLim',[-0.18,0.18],... % force y-axis to good limits
          'linFit','ON');        % show linear fit
  
      
end

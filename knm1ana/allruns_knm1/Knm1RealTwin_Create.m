function [Real , Twin] = Knm1RealTwin_Create(varargin)

p = inputParser;
p.addParameter('RealFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TwinFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('RunListFlag','KNM1');
p.parse(varargin{:});
RealFlag        = p.Results.RealFlag;
TwinFlag        = p.Results.TwinFlag;
RunListFlag     = p.Results.RunListFlag;

%% Option
options = {...
    'RunList',RunListFlag,...
    'exclDataStart',14,...
    'fixPar','1 5 6 7 8 9 10 11',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinT2',...
    'StackTolerance',1,...
    'NonPoissonScaleFactor',1.064,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]};

%% Definition of a Data MultiRunAnalysis
switch RealFlag
    case 'ON'
        Real=MultiRunAnalysis('DataType','Real',options{:});
    case 'OFF'
        if ~exist('Twin','var')
            Real = [];
        end
end

%% Definition of a TWIN MultiRunAnalysis
switch TwinFlag
    case 'ON'
        Twin=MultiRunAnalysis(...
            'DataType','Twin',options{:});
    case 'OFF'
        if ~exist('Twin','var')
            Twin = [];
        end
end

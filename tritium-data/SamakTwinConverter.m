function SamakTwinConverter(varargin)
% Convert Samak Monte Carlo data to hdf5-file format
% 1. Copy file with real data to TwinSamak directory (already in hdf5-format)
% 2. Replace counts in this file with counts from Samak twin file
% Lisa, January 2020

p = inputParser;
p.addParameter('DataSet','Knm2',@(x)ischar(x));
p.addParameter('RunNr',56341,@(x)isfloat(x));
p.addParameter('TwinLabel','_E018573.70eV_WGTSMACE_NIS1',@(x)ischar(x));
p.parse(varargin{:});
DataSet            = p.Results.DataSet;
RunNr              = p.Results.RunNr;
TwinLabel          = p.Results.TwinLabel;

% define which runsummary version to load
switch DataSet
    case 'FirstTritium.katrin'
        Version = '4a-fpd00'; % runn = 40531:41033;
    case 'Knm1'
        Version = 'RunSummary-Durable3a-fpd00';
    case 'Twin_Fitrium_Knm1'
        %version = 'Fitrium';
        Version = 'RunSummary-Durable2a-fpd00';
    case  'Twin_Kafit_Knm1'
        Version = 'MCRunSummary_Run00';
    case 'Knm2'
        Version = 'RunSummary-Prompt4b-fpd00';
end

%get the HD5 files copy them in the new folder and rename the files in the new folder
h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'/']; % path of data hdf5 format
h5name = [Version,num2str(RunNr),'.h5'];
    
% Create new directory for Samak twins in hdf5 format
twinh5path = [getenv('SamakPath'),'tritium-data/hdf5/',DataSet,'_TwinSamak/'];
MakeDir(twinh5path);
if contains(TwinLabel,'mNuSq')
    extraStr = 'mNu_';
else
    extraStr = '';
end
if contains(TwinLabel,'WGTSMACE_NIS1')
    extraStr = [extraStr,'WGTSMACE_NIS1_'];
elseif contains(TwinLabel,'MACE')
    extraStr = [extraStr,'MACE_'];
end

twinh5name = ['Twin_',extraStr,h5name];

% Copy hdf5 (real) data file to Samak hdf5 twin folder
copyfile([h5path,h5name],[twinh5path,twinh5name]);

% load Samak twins (mat format)
twinmatpath =[getenv('SamakPath'), 'tritium-data/mat/Twin',DataSet,'/'];
twinmatname = ['Twin',num2str(RunNr),TwinLabel,'.mat'];
TwinData = load([twinmatpath,twinmatname],'TBDIS');

% add -300eV point for KNM2
if strcmp(DataSet,'Knm2')
    TwinTBDIS = [zeros(1,148);TwinData.TBDIS];
end

% get correct ordering from RS
RealDataqU = h5read([h5path,h5name],'/RunSummary/RetardingEnergies/K35Readings');
[~, qUIndexSort] = sort(RealDataqU);               % Index from unsorted to sorted
qUIndexUnsort(qUIndexSort) = 1:numel(qUIndexSort); % Index from sorted to unsorted
TwinTBDIS = TwinTBDIS(qUIndexUnsort,:)';           % sort Twin TBDIS according to order in run-summary (unsorted)

% set efficiency correction to 1 -> NO corrections for twins
EffCorr   = ones(148,size(RealDataqU,1));

% Replace counts in hdf5 file with counts from Samak twins
h5write([twinh5path,twinh5name],'/RunSummary/Counts',TwinTBDIS);
h5write([twinh5path,twinh5name],'/RunSummary/EfficiencyCorrections',EffCorr);

fprintf('Convert Samak twin to hdf5 - save file to %s \n',[twinh5path,twinh5name]);
end




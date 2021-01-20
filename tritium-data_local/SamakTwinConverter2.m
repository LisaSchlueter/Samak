function SamakTwinConverter2(varargin)

% function reading Knm1 runs hd5 files copy them and replace counts with
% Samak Twin counts
%Lisa Schlueter - Thierry Lasserre -Joel Dai
% Last Modified: June 2019
%startup;

% Parser
p = inputParser;

p.addParameter('DataSet','Knm2',@(x)ischar(x));
p.addParameter('SaveConvertedFiles','ON',@(x)ismember(x,{'ON','OFF'})); %
p.addParameter('RunNr',56341,@(x)isfloat(x));

p.parse(varargin{:});

DataSet            = p.Results.DataSet;
SaveConvertedFiles = p.Results.SaveConvertedFiles;
RunNr              = p.Results.RunNr;
% 
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
%select the right files
h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'/'];     % path of data hdf5 format
h5name = [h5path,Version,num2str(RunNr),'.h5'];
    
% Create new directory for Samak twins in hdf5 format
twinh5path = [getenv('SamakPath'),'tritium-data/hdf5/TwinSamak_',DataSet,'/'];

%MakeDir(twinh5path);

% Copy hdf5 (real) data file to Samak hdf5 twin folder
%copyfile();


end
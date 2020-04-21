function HDF5readallruns(varargin)
p=inputParser;
p.addParameter('DataSet','Knm2',@(x)ischar(x));
p.addParameter('reConvert','ON',@(x)ismember(x,{'ON','OFF'})); % convert all files (again), which are in h5 folder even though they may be already converted
p.addParameter('h5runlist',[],@(x)isfloat(x));
p.parse(varargin{:});
DataSet = p.Results.DataSet;
reConvert = p.Results.reConvert;
h5runlist = p.Results.h5runlist;

notex = 1; ex = 2;

switch DataSet
    case 'FirstTritium.katrin'
        version = '4a-fpd00'; % runn = 40531:41033;
    case 'Knm1'
        version = 'RunSummary-Durable3a-fpd00';
    case 'Twin_Fitrium_Knm1'
        %version = 'Fitrium';
        version = 'RunSummary-Durable2a-fpd00';
    case  'Twin_Kafit_Knm1'
        version = 'MCRunSummary_Run00';
    case 'Knm2'
        version = 'RunSummary-Durable5d-fpd00';
end

if isempty(h5runlist)
% build list of runs in h5 folder
tmp = dir([getenv('SamakPath'),sprintf('/tritium-data/hdf5/%s',DataSet)]);
h5list = arrayfun(@(x) x.name,tmp,'UniformOutput',0);
h5list(~contains(h5list,'.h5'))=[];
h5list(~contains(h5list,version))=[];
if contains(DataSet,'Fitrium')
  h5list =  extractBefore(h5list,'_mc.h5');
else
h5list =  extractBefore(h5list,'.h5');
end
h5runlist = str2double(extractAfter(h5list,version));
end

% build list of runs in mat folder
tmp = dir([getenv('SamakPath'),sprintf('/tritium-data/mat/%s',DataSet)]);
matlist = arrayfun(@(x) strrep(x.name,version,''),tmp,'UniformOutput',0);
matlist(~contains(matlist,'.mat'))=[];
matrunlist = strrep(matlist,'.mat','');
matrunlist = cell2mat(cellfun(@(x) str2double(x),matrunlist,'UniformOutput',0));
matrunlist = matrunlist(~isnan(matrunlist));

if strcmp(reConvert,'OFF')% select only runs, which don't exist in mat folder
    runn =h5runlist(~ismember(h5runlist,matrunlist));
else % convert all (again)
    runn = h5runlist;
end

if isempty(runn)
    fprintf(2,'no new runs to convert \n');
    return
end
% start conversion

for ii = 1:length(runn)
    try
        switch DataSet
            case {'Knm1','FirstTritium.katrin','Knm2'}
                HDF5Reader('RunNr',runn(ii),'version',version);
            case 'Twin_Fitrium_Knm1'
                HDF5Reader('RunNr',runn(ii),'Fitter','Fitrium','version',version);
            case 'Twin_Kafit_Knm1'
                 HDF5Reader('RunNr',runn(ii),'Fitter','Kafit','version',version);
        end 
        existing(ex) = runn(ii);
        fprintf('run %0.0d converted to mat file \n',runn(ii));
        ex = ex + 1;
    catch ErrorMessage
        if contains(ErrorMessage.identifier,'fileOpenErr')
           fprintf('File doesnt exist: Run %0.0d \n',runn(ii))
        else
             fprintf('Problem with Run %0.0d: \n',runn(ii))
             ErrorMessage.message
        end
       notexisting(notex) = runn(ii);
        notex = notex + 1;
    end
end
end

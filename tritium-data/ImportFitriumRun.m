function  ImportFitriumRun(varargin)

startup;
% Parser 
p = inputParser;
p.addParameter('DataSet','Knm1',@(x)ischar(x));
p.addParameter('Username','jdai',@(x)ischar(x));
p.parse(varargin{:});
DataSet            = p.Results.DataSet;
Username           = p.Results.Username;

twinpath = [getenv('SamakPath'), 'tritium-data/mat/Twin',DataSet,'/'];
folderpath = [getenv('SamakPath'), 'tritium-data/hdf5/Twin_Fitrium_',DataSet,'/'];
NbFiles = length(dir(twinpath));
TwinObj = struct2cell(dir(twinpath));
twinnamelist = string(TwinObj(1,:));
progressbar('Getting Twin ....');
for i=1:NbFiles
    progressbar(i/NbFiles)
    if contains(twinnamelist(i),'Twin')
        twinrunnumber = extractBefore(twinnamelist(i),'.mat');
        twinrunnumber = char(extractAfter(twinrunnumber,'Twin'));
        system(['scp ',Username,'@pcltr-01.mpp.mpg.de:../karlch/Fitting/KNM1/MC/MCDataNoFluct/RunSummary-Durable2a-fpd00',twinrunnumber,'_mc.h5 Fitrium',twinrunnumber,'.h5']);
        system(['mv ~/Samak2.0/Fitrium',twinrunnumber,'.h5 ',folderpath ]);
    end
end
end



function [h5mat] = HDF5Reader(varargin)
% Function reading .ktf KATRIN Run Summary Files
% Convert into HDF5 Files
% Lisa Schlueter - Thierry Lasserre
% Last Modified: June 2019

% Parser 
p = inputParser;
% Main inputs (necessary)
p.addParameter('Version','3c-fpd00',@(x) ischar(x)); %Version number of the Run Summaries
p.addParameter('RunNr',40667,@(x)all(isfloat(x)));
p.addParameter('Fitter','Samak',@(x)ismember(x,{'Samak','Fitrium','Kafit'}));
p.addParameter('TimeBias',1,@(x)isfloat(x)); % increase measurement time by natural number (for Fitritum twins)
p.addParameter('ExtraLabel','',@(x)ischar(x)); % extra label for Fitrium twins
%p.addParameter('SamakPath','../../');
% optional inputs
p.addParameter('saveTD_DataBank','OFF',@(x)ismember(x,{'ON','OFF'})); %Save created TD per Run (not neccessary, also included in mat file)
p.parse(varargin{:});

Version         = p.Results.Version;
RunNr           = p.Results.RunNr;
saveTD_DataBank = p.Results.saveTD_DataBank;
Fitter          = p.Results.Fitter;
TimeBias        = p.Results.TimeBias;
ExtraLabel      = p.Results.ExtraLabel;

% define data set label
DataSet = GetDataSet(RunNr);
switch Fitter
    case 'Samak'
        h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'/'];
        h5name = [Version,num2str(RunNr),'.h5'];
    case 'Fitrium'
        h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'_TwinFitrium/'];
        h5name = [ExtraLabel,Version,num2str(RunNr),'.h5'];
        samakh5name = [Version,num2str(RunNr),'.h5'];
        samakh5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'/'];
    case 'Kafit'
        h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'_TwinKaFit/'];
        h5name = [ExtraLabel,Version,num2str(RunNr),'.h5'];
        samakh5name = [Version,num2str(RunNr),'.h5'];
        samakh5path = [getenv('SamakPath'), 'tritium-data/hdf5/',DataSet,'/'];
end
%% Import Data Counts
TBDISuncorr = h5read([h5path,h5name],'/RunSummary/Counts')'; % for compute TBDISE
try
    EffCorr     = double(h5read([h5path,h5name],'/RunSummary/EfficiencyCorrections')');
catch
    EffCorr     = TBDISuncorr./TBDISuncorr;
end
EffCorr(EffCorr==0) = 1;
TBDIS       = double(TBDISuncorr)./double(EffCorr).*TimeBias;
TBDISE      = sqrt(TBDIS./EffCorr); % statstical uncertainty

if strcmp(DataSet,'Knm2') && strcmp(Fitter,'Samak')
    if contains(Version,'Durable')
        TBDIS1133 = h5read([h5path,h5name],'/RunSummary/Counts1133')'; % for compute TBDISE
        TBDIS2232 = h5read([h5path,h5name],'/RunSummary/Counts2232')'; % for compute TBDISE
        TBDIS14keV = zeros(size(TBDIS));
    else
        TBDIS14keV = h5read([h5path,h5name],'/RunSummary/CountsKNM1')'; % for compute TBDISE
        TBDIS1133 = zeros(size(TBDIS));
        TBDIS2232 = zeros(size(TBDIS));
    end
else
    TBDIS14keV = zeros(size(TBDIS));
    TBDIS1133 = zeros(size(TBDIS));
    TBDIS2232 = zeros(size(TBDIS));
end
%% Import Retarding Potentials and Measurement Times
if strcmp(DataSet,'FirstTritium.katrin') && str2double(Version(1))<4
        qU = h5read([h5path,h5name],'/RunSummary/RetardingEnergy');
else
        K35       = h5read([h5path,h5name],'/RunSummary/RetardingEnergies/K35Readings'); %VnqU
        qUOffsets = h5read([h5path,h5name],'/RunSummary/RetardingEnergies/RetardingEnergyOffsets'); %V148
        tmpA      = repmat(K35',numel(qUOffsets),1);
        tmpB      = repmat(qUOffsets,1,numel(K35));
        qU        = tmpA+tmpB;
end

qU = qU';
TimeperSubRunperPixel = h5read([h5path,h5name],'/RunSummary/LiveTime')';
qUfrac = TimeperSubRunperPixel./sum(TimeperSubRunperPixel,1);
TimeSec = sum(TimeperSubRunperPixel,1)'.*TimeBias;
%% Import pixelwise magnetic field corrections (Ba)
if strcmp(GetDataSet(RunNr),'FirstTritium.katrin') && str2double(Version(1))<4
    MACE_Ba_T = 6.*ones(148,1).*1e-04;
    StartTimeStamp = 0;
else% RunNr >= 51021
    MACE_Ba_T = h5read([h5path,h5name],'/RunSummary/Transmission/MagneticFieldAna');
    switch Fitter
        case 'Samak'
            StartTimeStamp = h5read([h5path,h5name],'/RunSummary/Metadata/StartTimestamp');
        case {'Fitrium','Kafit'}
            StartTimeStamp = h5read([samakh5path,samakh5name],'/RunSummary/Metadata/StartTimestamp');
    end
end

%% Pinch magnetic field (Bmax)
VersionTmp = extractAfter(Version,'Prompt');
if isempty(VersionTmp)
 VersionTmp = extractAfter(Version,'Durable');   
end
if ismember(DataSet,{'FirstTritium.katrin','Knm1'})
    MACE_Bmax_T = 4.23.*ones(148,1); % not available in RS
else%if strcmp(DataSet,'Knm2') && str2double(VersionTmp(1))>=4
    MACE_Bmax_T = h5read([h5path,h5name],'/RunSummary/Transmission/MagneticFieldMax');
end
%% Import Slow Control Parameter:  Column Density and cross section

WGTS_CD_MolPerM2 = h5read([h5path,h5name],'/RunSummary/Source/ColumnDensityMean');
WGTS_CD_MolPerCm2 = WGTS_CD_MolPerM2*1e-4;

if ~strcmp(GetDataSet(RunNr),'FirstTritium.katrin') && ~ismember(Fitter,{'Fitrium','Kafit'})
    WGTS_CD_MolPerM2_Sigma =  h5read([h5path,h5name],'/RunSummary/Source/ColumnDensityCrossSectionMean');
    ISXsection = WGTS_CD_MolPerM2_Sigma./WGTS_CD_MolPerM2;
else
    ISXsection = 0;
end

% Column Density per subrun (not always available)
try
    if strcmp(Fitter,'Kafit')
         WGTS_CD_MolPerCm2_SubRun = repmat(WGTS_CD_MolPerCm2,size(qU,1),1);
    else
    WGTS_CD_MolPerCm2_SubRun = h5read([h5path,h5name],'/RunSummary/Source/ColumnDensity')*1e-4;   
    end
catch
    WGTS_CD_MolPerCm2_SubRun = repmat(WGTS_CD_MolPerCm2,numel(qU),1);
end

%% WGTS LARA Data: Isotopologue concentration
WGTS_MolFrac_TT = h5read([h5path,h5name],'/RunSummary/Source/LARA/T2FractionMean');
WGTS_MolFrac_DT = h5read([h5path,h5name],'/RunSummary/Source/LARA/DTFractionMean');
WGTS_MolFrac_HT = h5read([h5path,h5name],'/RunSummary/Source/LARA/HTFractionMean');

try
    % LARA Data per subrun
    if strcmp(Fitter,'Kafit')
        WGTS_MolFrac_TT_SubRun = repmat(WGTS_MolFrac_TT,size(qU,1),1);
        WGTS_MolFrac_DT_SubRun = repmat(WGTS_MolFrac_DT,size(qU,1),1);
        WGTS_MolFrac_HT_SubRun = repmat(WGTS_MolFrac_HT,size(qU,1),1);
    else
    WGTS_MolFrac_TT_SubRun = h5read([h5path,h5name],'/RunSummary/Source/LARA/T2Fraction');
    WGTS_MolFrac_DT_SubRun = h5read([h5path,h5name],'/RunSummary/Source/LARA/DTFraction');
    WGTS_MolFrac_HT_SubRun = h5read([h5path,h5name],'/RunSummary/Source/LARA/HTFraction');
    end
catch
    WGTS_MolFrac_TT_SubRun=repmat(WGTS_MolFrac_TT,numel(qU),1);
    WGTS_MolFrac_DT_SubRun=repmat(WGTS_MolFrac_DT,numel(qU),1);
    WGTS_MolFrac_HT_SubRun=repmat(WGTS_MolFrac_HT,numel(qU),1);
end
try
    % LARA Data error per subrun
    WGTS_MolFrac_TT_SubRun_error = h5readatt(h5name,'/RunSummary/Source/LARA/T2Fraction','Error');
    WGTS_MolFrac_DT_SubRun_error = h5readatt(h5name,'/RunSummary/Source/LARA/DTFraction','Error');
    WGTS_MolFrac_HT_SubRun_error = h5readatt(h5name,'/RunSummary/Source/LARA/HTFraction','Error');
catch
    WGTS_MolFrac_TT_SubRun_error=zeros(size(qU,1),1);
    WGTS_MolFrac_DT_SubRun_error=zeros(size(qU,1),1);
    WGTS_MolFrac_HT_SubRun_error=zeros(size(qU,1),1);
end

WGTS_MolFrac_TT_SubRun_error_mean = sqrt(sum(WGTS_MolFrac_TT_SubRun_error.^2))./length(WGTS_MolFrac_TT_SubRun_error);
WGTS_MolFrac_DT_SubRun_error_mean = sqrt(sum(WGTS_MolFrac_DT_SubRun_error.^2))./length(WGTS_MolFrac_DT_SubRun_error);
WGTS_MolFrac_HT_SubRun_error_mean = sqrt(sum(WGTS_MolFrac_HT_SubRun_error.^2))./length(WGTS_MolFrac_HT_SubRun_error);

%% import Rear wall information (available since KNM2 - version Prompt 4b)

if strcmp(DataSet,'Knm2') && str2double(VersionTmp(1))>=4 && strcmp(Fitter,'Samak')
    RW_BiasVoltage = h5read([h5path,h5name],'/RunSummary/RearWall/BiasVoltageMean');
else
    RW_BiasVoltage = 0;
end
%% special treatment for run with 39 subruns (background subrun is missing)
if RunNr==51861
    WGTS_CD_MolPerCm2_SubRun     = [0;WGTS_CD_MolPerCm2_SubRun];
    WGTS_MolFrac_TT_SubRun       = [0;WGTS_MolFrac_TT_SubRun];
    WGTS_MolFrac_HT_SubRun       = [0;WGTS_MolFrac_HT_SubRun];
    WGTS_MolFrac_DT_SubRun       = [0;WGTS_MolFrac_DT_SubRun];
    WGTS_MolFrac_TT_SubRun_error = [0;WGTS_MolFrac_TT_SubRun_error];
    WGTS_MolFrac_DT_SubRun_error = [0;WGTS_MolFrac_DT_SubRun_error];
    WGTS_MolFrac_HT_SubRun_error = [0;WGTS_MolFrac_HT_SubRun_error];
end
%% rearrage subrun values if neccesary
% subruns arrays are always arranged as upward scans
if  (any(diff(qU(:,1)) < 0))
    [qU,SortedIndex] = sort(qU);
    TBDIS      = TBDIS(SortedIndex(:,1),:);
    qUfrac     = qUfrac(SortedIndex(:,1),:);
    EffCorr    = EffCorr(SortedIndex(:,1),:);
    TBDIS14keV = TBDIS14keV(SortedIndex(:,1),:);
    TBDIS1133  = TBDIS1133(SortedIndex(:,1),:);
    TBDIS2232  = TBDIS2232(SortedIndex(:,1),:);
   
    TimeperSubRunperPixel = TimeperSubRunperPixel(SortedIndex(:,1),:);
    if ~strcmp(Fitter,'Kafit')
        WGTS_CD_MolPerCm2_SubRun = WGTS_CD_MolPerCm2_SubRun(SortedIndex(:,1));
        WGTS_MolFrac_TT_SubRun = WGTS_MolFrac_TT_SubRun(SortedIndex(:,1));
        WGTS_MolFrac_DT_SubRun = WGTS_MolFrac_DT_SubRun(SortedIndex(:,1));
        WGTS_MolFrac_HT_SubRun = WGTS_MolFrac_HT_SubRun(SortedIndex(:,1));
        WGTS_MolFrac_TT_SubRun_error = WGTS_MolFrac_TT_SubRun_error(SortedIndex(:,1));
        WGTS_MolFrac_DT_SubRun_error = WGTS_MolFrac_DT_SubRun_error(SortedIndex(:,1));
        WGTS_MolFrac_HT_SubRun_error = WGTS_MolFrac_HT_SubRun_error(SortedIndex(:,1));
    end
end

%% Temperature
if strcmp(Fitter,'Samak')
    WGTS_Temp_K = h5read(h5name,'/RunSummary/Source/TemperatureMean');
end
%% Special Treatment for some Data Points
% First Tritium (FT): Erase two first points from the two first runs

switch DataSet
    case 'FirstTritium.katrin'
        if (RunNr == 40257) || (RunNr == 40258)
            qU = qU(3:28,:);
            qUfrac = qUfrac(3:28,:);
            TBDIS = TBDIS(3:28,:);
            TimeperSubRunperPixel = TimeperSubRunperPixel(3:28,:);
        end
        
        qU_RM                       = 0;
        qUfrac_RM                   = 0;
        TBDIS_RM                    = 0;
        TimeperSubRunperPixel_RM    = 0;
    case {'Knm2','Knm1'}%{'Knm1','Twin_Fitrium_Knm1','Twin_Kafit_Knm1','Knm2'}
        qU_RM                       = qU(1,:)';
        qUfrac_RM                   = qUfrac(1,:)';
        TBDIS_RM                    = TBDIS(1,:)';
        EffCorr_RM                  =  EffCorr(1,:)';
        TBDIS14keV_RM               = TBDIS14keV(1,:)';
        TBDIS1133_RM                = TBDIS1133(1,:)';
        TBDIS2232_RM                = TBDIS2232(1,:)';
        TimeperSubRunperPixel_RM    = TimeperSubRunperPixel(1,:)';
        
        % delete rate monitor point: not used for nu-mass analysis
        qUstart = 2;
        qU                           = qU(qUstart:end,:);
        qUfrac                       = qUfrac(qUstart:end,:);
        TBDIS                        = TBDIS(qUstart:end,:);
        TBDIS14keV                   = TBDIS14keV(qUstart:end,:);
        TBDIS1133                    = TBDIS1133(qUstart:end,:);
        TBDIS2232                    = TBDIS2232(qUstart:end,:);
        TimeperSubRunperPixel        = TimeperSubRunperPixel(qUstart:end,:);
        EffCorr                      = EffCorr(qUstart:end,:);
        WGTS_CD_MolPerCm2_SubRun     = WGTS_CD_MolPerCm2_SubRun(qUstart:end);
        WGTS_MolFrac_TT_SubRun       = WGTS_MolFrac_TT_SubRun(qUstart:end);
        WGTS_MolFrac_DT_SubRun       = WGTS_MolFrac_DT_SubRun(qUstart:end);
        WGTS_MolFrac_HT_SubRun       = WGTS_MolFrac_HT_SubRun(qUstart:end);     
        WGTS_MolFrac_TT_SubRun_error = WGTS_MolFrac_TT_SubRun_error(qUstart:end);
        WGTS_MolFrac_DT_SubRun_error = WGTS_MolFrac_DT_SubRun_error(qUstart:end);
        WGTS_MolFrac_HT_SubRun_error = WGTS_MolFrac_HT_SubRun_error(qUstart:end);     
end


%% Save slow control data and counts PIXELWISE
save_path = strrep(h5path,'hdf5','mat');

if ~exist(save_path,'dir')
    system(['mkdir ',save_path])
end

switch Fitter
    case 'Samak'
        matFileName = [num2str(RunNr),'.mat'];
    case 'Fitrium'
        matFileName = ['TwinFitrium_',strrep(ExtraLabel,'mc_',''),num2str(RunNr),'.mat'];
    case 'Kafit'
        matFileName = ['TwinKafit_',strrep(ExtraLabel,'Kafit_',''),num2str(RunNr),'.mat'];
end

if TimeBias>1
    matFileName = strrep(matFileName,'.mat',sprintf('_%.0ftimes.mat',TimeBias));
end
MakeDir(save_path);
savename = [save_path,matFileName];

save(savename,'StartTimeStamp','TBDIS','TBDISE','qU','qUfrac','EffCorr',...
    'TimeSec','TimeperSubRunperPixel',...
    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
    'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun',...
    'WGTS_MolFrac_TT_SubRun_error','WGTS_MolFrac_HT_SubRun_error',...
    'WGTS_MolFrac_DT_SubRun_error','WGTS_MolFrac_TT_SubRun_error_mean',...
    'WGTS_MolFrac_DT_SubRun_error_mean','WGTS_MolFrac_HT_SubRun_error_mean',...
    'MACE_Ba_T','StartTimeStamp','ISXsection',...
    'qU_RM','qUfrac_RM','TBDIS_RM','EffCorr_RM',...
    'RW_BiasVoltage',...
    'MACE_Bmax_T',...
    'matFileName',...
    'TimeBias',...
    'TBDIS14keV','TBDIS14keV_RM',...
    'TBDIS1133','TBDIS1133_RM',...
    'TBDIS2232','TBDIS2232_RM',...
    'WGTS_Temp_K',...
    '-v7.3','-nocompression');

h5mat.mpix = struct('StartTimeStamp',StartTimeStamp,'TBDIS',TBDIS,'TBDISE',TBDISE,'qU',qU,'TimeperSubRunperPixel',TimeperSubRunperPixel,...
    'TimeSec',TimeSec,'qUfrac',qUfrac,...
    'EffCorr',EffCorr,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_CD_MolPerCm2_SubRun', WGTS_CD_MolPerCm2_SubRun,...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...
    'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...
    'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
    'WGTS_MolFrac_TT_SubRun_error',WGTS_MolFrac_TT_SubRun_error,...
    'WGTS_MolFrac_DT_SubRun_error',WGTS_MolFrac_DT_SubRun_error,...
    'WGTS_MolFrac_HT_SubRun_error',WGTS_MolFrac_HT_SubRun_error,...
    'WGTS_MolFrac_TT_SubRun_error_mean',WGTS_MolFrac_TT_SubRun_error_mean,...
    'WGTS_MolFrac_DT_SubRun_error_mean',WGTS_MolFrac_DT_SubRun_error_mean,...
    'WGTS_MolFrac_HT_SubRun_error_mean',WGTS_MolFrac_HT_SubRun_error_mean,...
    'MACE_Ba_T',MACE_Ba_T,...
    'qU200',qU_RM,'qUfrac200',qUfrac_RM,'TBDIS200',TBDIS_RM,...
    'ISXsection',ISXsection,...
    'RW_BiasVoltage',RW_BiasVoltage,...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeBias',TimeBias,....
    'TBDIS14keV',TBDIS14keV,...
    'TBDIS14keV_RM',TBDIS14keV_RM,...
    'TBDIS1133',TBDIS1133,...
    'TBDIS1133_RM',TBDIS1133_RM,...     
    'TBDIS2232',TBDIS2232,...
    'TBDIS2232_RM',TBDIS2232_RM,...
    'WGTS_Temp',WGTS_Temp_K);

% TD (optional)
switch saveTD_DataBank
    case 'ON'
        RunTime = TimeSec;
        TD = ['Run',num2str(RunNr),'mpix'];
            save([SamakPath '/simulation/katrinsetup/TD_DataBank/Run',...
            num2str(RunNr),'mpix','.mat'],...
            'qU','qUfrac','RunTime','TD','-v7.3','-nocompression');
        
        %Stacked TD
        qU = mean(qU,2);
        
        TD = ['Run',num2str(RunNr)];
                 save([SamakPath '/simulation/katrinsetup/TD_DataBank/Run',...
            num2str(RunNr),'.mat'],...
            'qU','qUfrac','RunTime','TD','-v7.3','-nocompression');
end
end


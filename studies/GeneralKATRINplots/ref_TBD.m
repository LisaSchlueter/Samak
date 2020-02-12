function TritiumObject = ref_RunAnalysis_Knm1(DataType,D,RunNr_or_StackFileName,TwinName, PixList,RingList,varargin)%,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION FILE - Very First Tritium
% For all kind of analysis
% Single Runs:
% - Stacked pixels (Uniform Mode)
% Stacked Runs:
% - Stacked pixels (Uniform Mode)
% - Ring
% - Pixels

% It can recieve as input in RunNr_StackFileName either a run number or
% the file name of a stacked runs file.
%
% T. Lasserre 2018
% L. Schlueter
% P. I. Morales Guzman
% Thierry Lasserre
% Last update 15/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------- %
% Label for TD & TBD Simulation Object:
switch DataType
    % BASED ON REAL DATA RUNS
    case 'Real'
        TDMode = 'DataTBD'; 
        RunPrefix='';
    case 'Twin'
        TDMode = 'DataTBD'; 
        RunPrefix = 'Twin';
    case 'Fake'
        %        RunPrefix = 'Fake1Run';
        %        Fake1RunSimObj=load('Fake1ex2RunObject.mat');
        %        RunPrefix = 'Fake5yRun';
        %        Fake1RunSimObj=load('Fake5yex2RunObject.mat');
        TDMode = 'Read';
        RunPrefix              = 'FakeKITNuScanRun';
        Fake1RunSimObj         = load('FakeKITNuScanex2RunObject.mat');
end
if isfloat(RunNr_or_StackFileName)
    RunSummaryFileName = [RunPrefix,num2str(RunNr_or_StackFileName),TwinName,'.mat'];
elseif ischar(RunNr_or_StackFileName)
    RunSummaryFileName = [RunNr_or_StackFileName,'.mat'];
end
% FAKE
%   D = load(RunSummaryFileName); % for SingleRunObj: take Data from RunSummary
% ---------------------------------------------------------------------- %

% ---------------------------------------------------------------------- %
% Parser Parameters from RunSummary
p = inputParser;
p.addParameter('TimeSec',D.TimeSec,@(x)isfloat(x));
switch DataType
    case {'Real','Twin'}
        %p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2/0.978,@(x)isfloat(x) && x>0);%4.560143569452563e+17 % WARNING First Tritium!!!!!
%        p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2*1.046,@(x)isfloat(x) && x>0); %4.560143569452563e+17 % WARNING KNM1 !!!!!
        p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2,@(x)isfloat(x) && x>0); %4.560143569452563e+17 % WARNING!!!!!
    case 'Fake'
        p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2,@(x)isfloat(x) && x>0);% 4.560143569452563e+17 % WARNING!!!!!
end
p.addParameter('WGTS_CD_MolPerCm2_SubRun',D.WGTS_CD_MolPerCm2_SubRun);
p.addParameter('WGTS_MolFrac_TT',D.WGTS_MolFrac_TT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT_SubRun',D.WGTS_MolFrac_TT_SubRun);
p.addParameter('WGTS_MolFrac_DT',D.WGTS_MolFrac_DT,@(x)isfloat(x) && x>0); 
p.addParameter('WGTS_MolFrac_DT_SubRun',D.WGTS_MolFrac_DT_SubRun);
p.addParameter('WGTS_MolFrac_HT',D.WGTS_MolFrac_HT,@(x)isfloat(x)&& x>0);
p.addParameter('WGTS_MolFrac_HT_SubRun',D.WGTS_MolFrac_HT_SubRun);
p.addParameter('qU',D.qU,@(x)isfloat(x) && all(x>0)); % if empty, read from data files
p.addParameter('qUfrac',D.qUfrac,@(x)isfloat(x));
            
% MTD
SplitStringCell = regexp(RunSummaryFileName, '/', 'split');
p.addParameter('TD',strrep(SplitStringCell{end},'.mat',''));

switch DataType 
    case {'Real','Twin'}
        % Binning
        p.addParameter('nTeBinningFactor',100,@(x)isfloat(x) && x>0);
        % KATRIN GENERAL SETTINGS
        p.addParameter('WGTS_B_T',2.52,@(x)isfloat(x) && x>0);
        p.addParameter('MACE_Bmax_T',4.23,@(x)isfloat(x) && x>0);
        p.addParameter('MACE_Ba_T',D.MACE_Ba_T,@(x)isfloat(x) && x>0);
        p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));
        % TD Reading Mode
        p.addParameter('mnuSq_i',0,@(x)isfloat(x));
        % ---------------------------------------------------------------------- %
        % FT
        % p.addParameter('Q_i',18573.7,@(x)isfloat(x)); % effective Endpoint
        % ---------------------------------------------------------------------- %
        % KITNUSCAN
        p.addParameter('Q_i',18573.7,@(x)isfloat(x)); % effective Endpoint
        % ---------------------------------------------------------------------- %
        % WGTS: Flags for FSD: T-T / D-T / H-T  
        p.addParameter('DTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','BlindingKNM1'})); % KITNUSCAN - Blinded FSD
        p.addParameter('HTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','SAENZ','BlindingKNM1'}));
        p.addParameter('TTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','SAENZ','BlindingKNM1'}));
        % ---------------------------------------------------------------------- %
        p.addParameter('ELossFlag','KatrinD2',@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','KatrinD2'}));
        p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));
        % FPD
        p.addParameter('FPD_MeanEff',0.9,@(x)isfloat(x) && x>0);
        p.addParameter('FPD_ROIlow',14,@(x)isfloat(x) && x>0);
        p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
        p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
        p.addParameter('BKG_RateAllFPDSec','');
        p.addParameter('BKG_RatePixelSec','');
        p.addParameter('BKG_RateRingSec','');
        p.addParameter('FPD_Segmentation','OFF',@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));
        % Transmission Function
        p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));
        p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('UseParallelRF','ON',@(x)ismember(x,{'OFF','ON'}));
    
    case 'Fake'
        % Binning
        p.addParameter('nTeBinningFactor',Fake1RunSimObj.SimFakeObj_i.nTeBinningFactor,@(x)isfloat(x) && x>0);
        % KATRIN GENERAL SETTINGS
        p.addParameter('WGTS_B_T',Fake1RunSimObj.SimFakeObj_i.WGTS_B_T,@(x)isfloat(x) && x>0);
        p.addParameter('MACE_Bmax_T',Fake1RunSimObj.SimFakeObj_i.MACE_Bmax_T,@(x)isfloat(x) && x>0);
        p.addParameter('MACE_Ba_T',Fake1RunSimObj.SimFakeObj_i.MACE_Ba_T,@(x)isfloat(x) && x>0);
        p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));
        % TD Reading Mode
        p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));
        p.addParameter('mnuSq_i',Fake1RunSimObj.SimFakeObj_i.mnuSq_i,@(x)isfloat(x));
        p.addParameter('Q_i',Fake1RunSimObj.SimFakeObj_i.Q_i,@(x)isfloat(x)); % effective Endpoint
        % WGTS: Flags for FSD: T-T / D-T / H-T
        % ---------------------------------------------------------------------- %
        % FT
        % p.addParameter('DTFSD',Fake1RunSimObj.SimFakeObj_i.DTFSD,@(x)ismember(x,{'OFF','DOSS','BlindingKNM1'}));
        % p.addParameter('HTFSD',Fake1RunSimObj.SimFakeObj_i.HTFSD,@(x)ismember(x,{'OFF','SAENZ','BlindingKNM1'}));
        % p.addParameter('TTFSD',Fake1RunSimObj.SimFakeObj_i.TTFSD,@(x)ismember(x,{'OFF','DOSS','SAENZ','BlindingKNM1'}));
        % ---------------------------------------------------------------------- %
        % TEST BLINDING FSD
        p.addParameter('DTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','BlindingKNM1'}));
        p.addParameter('HTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','SAENZ','BlindingKNM1'}));
        p.addParameter('TTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','SAENZ','BlindingKNM1'}));
        % ---------------------------------------------------------------------- %
        p.addParameter('ELossFlag',Fake1RunSimObj.SimFakeObj_i.ELossFlag,@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','KatrinD2'}));
        p.addParameter('DopplerEffectFlag',Fake1RunSimObj.SimFakeObj_i.DopplerEffectFlag,@(x)ismember(x,{'OFF','numConv','matConv'}));
        % FPD
        p.addParameter('FPD_MeanEff',Fake1RunSimObj.SimFakeObj_i.FPD_MeanEff,@(x)isfloat(x) && x>0);
        p.addParameter('FPD_ROIlow',Fake1RunSimObj.SimFakeObj_i.FPD_ROIlow,@(x)isfloat(x) && x>0);
        p.addParameter('FPD_ROIEff',Fake1RunSimObj.SimFakeObj_i.FPD_ROIEff,@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('FPD_PileUpEff',Fake1RunSimObj.SimFakeObj_i.FPD_PileUpEff,@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('BKG_Flag',Fake1RunSimObj.SimFakeObj_i.BKG_Flag,@(x)ismember(x,{'ON','OFF','XmasData'}));
        p.addParameter('BKG_Type',Fake1RunSimObj.SimFakeObj_i.BKG_Type,@(x)ismember(x,{'FLAT','SLOPE'}));
        p.addParameter('BKG_RateAllFPDSec',Fake1RunSimObj.SimFakeObj_i.BKG_RateAllFPDSec);
        p.addParameter('BKG_RateRingSec','');
        p.addParameter('FPD_Segmentation',Fake1RunSimObj.SimFakeObj_i.FPD_Segmentation,@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));
        % Transmission Function
        p.addParameter('KTFFlag',Fake1RunSimObj.SimFakeObj_i.KTFFlag,@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));
        p.addParameter('recomputeRF',Fake1RunSimObj.SimFakeObj_i.recomputeRF,@(x)ismember(x,{'ON','OFF'}));
        %p.addParameter('recomputeRF','ON',@(x)ismember(x,{'ON','OFF'}));
        p.addParameter('UseParallelRF',Fake1RunSimObj.SimFakeObj_i.UseParallelRF,@(x)ismember(x,{'OFF','ON'}));

end

p.parse(varargin{:});

%RunSummary Parameter
TimeSec                  = p.Results.TimeSec;
WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_CD_MolPerCm2_SubRun = p.Results.WGTS_CD_MolPerCm2_SubRun;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
WGTS_MolFrac_HT_SubRun   = p.Results.WGTS_MolFrac_HT_SubRun;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
WGTS_MolFrac_DT_SubRun   = p.Results.WGTS_MolFrac_DT_SubRun;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;
WGTS_MolFrac_TT_SubRun   = p.Results.WGTS_MolFrac_TT_SubRun;
TD                       = p.Results.TD;
qU                       = p.Results.qU;
qUfrac                   = p.Results.qUfrac;

% Parameters: Calculation precision
nTeBinningFactor         = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
WGTS_B_T                 = p.Results.WGTS_B_T;
MACE_Bmax_T              = p.Results.MACE_Bmax_T;
MACE_Ba_T                = p.Results.MACE_Ba_T; 
ISCS                     = p.Results.ISCS;
ELossFlag                = p.Results.ELossFlag;

FPD_MeanEff              = p.Results.FPD_MeanEff;
FPD_ROIlow               = p.Results.FPD_ROIlow;
FPD_ROIEff               = p.Results.FPD_ROIEff;
FPD_PileUpEff            = p.Results.FPD_PileUpEff;
BKG_Flag                 = p.Results.BKG_Flag;
BKG_Type                 = p.Results.BKG_Type;
BKG_RateAllFPDSec        = p.Results.BKG_RateAllFPDSec;
BKG_RateRingSec          = p.Results.BKG_RateRingSec;
BKG_RatePixelSec         = p.Results.BKG_RatePixelSec;
mnuSq_i                  = p.Results.mnuSq_i;
Q_i                      = p.Results.Q_i;
%FPD segmentation
FPD_Segmentation         = p.Results.FPD_Segmentation;

% Doppler Effect
DopplerEffectFlag        = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD                    = p.Results.DTFSD;
HTFSD                    = p.Results.HTFSD;
TTFSD                    = p.Results.TTFSD;

% Transmission Function 
KTFFlag                  = p.Results.KTFFlag;
recomputeRF              = p.Results.recomputeRF;
UseParallelRF            = p.Results.UseParallelRF;
% End Parser Parameters from RunSummary
% ---------------------------------------------------------------------- %

% ---------------------------------------------------------------------- %
% Create Tritium spectrum object 
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i,...
    'Q_i',Q_i,...
    'qU',qU,...
    'qUfrac',qUfrac};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...    
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
    'WGTS_MolFracRelErr_TT',std(WGTS_MolFrac_TT_SubRun)/WGTS_MolFrac_TT,...
    'WGTS_MolFracRelErr_DT',std(WGTS_MolFrac_DT_SubRun)/WGTS_MolFrac_DT,...
    'WGTS_MolFracRelErr_HT',std(WGTS_MolFrac_HT_SubRun)/WGTS_MolFrac_HT,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.5,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_CD_MolPerCm2_SubRun',WGTS_CD_MolPerCm2_SubRun,...    
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30,...
    'ISCS',ISCS,...
    'ELossFlag',ELossFlag,...
    'recomputeRF',recomputeRF,...
    'UseParallelRF',UseParallelRF,...
    'NIS',7}; %Thierry WARNING 13/2/2019 

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T,...
    'MACE_Ba_Setting','Data'}; %do not use pixelmap correction from txt file, use values from runsummary

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};


%% Get Background: Whole FPD (excluding some pixels), per Ring, per Pixel
[PixList,RingPixList] = Ring2Pixel(RingList,PixList);
if isempty(BKG_RateAllFPDSec)
    BKG_RateAllFPDSec = numel(PixList)/148*GetBackground('MACE_Ba_T',mean(MACE_Ba_T),...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
end
if isempty(BKG_RateRingSec) && strcmp(FPD_Segmentation,'RING')
    nPixperRing = cellfun(@(x) numel(x),RingPixList); nPixperRing(~ismember(1:13,RingList)) = [];
    BKG_RateRingSec = (nPixperRing./148)'.*GetBackground('MACE_Ba_T',MACE_Ba_T,...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
    BKG_RateRingSec(isnan(BKG_RateRingSec))=0; BKG_RateRingSec(isinf(BKG_RateRingSec))=0;
elseif isempty(BKG_RateRingSec) && ~strcmp(FPD_Segmentation,'RING')
    BKG_RateRingSec = BKG_RateAllFPDSec./numel(RingList).*ones(numel(RingList),1);
end

if isempty(BKG_RatePixelSec) && ismember(FPD_Segmentation,{'SINGLEPIXEL','MULTIPIXEL'})
    BKG_RatePixelSec = 1/148.*GetBackground('MACE_Ba_T',MACE_Ba_T,...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
    BKG_RatePixelSec(isnan(BKG_RatePixelSec))=0;  BKG_RatePixelSec(isinf(BKG_RatePixelSec))=0;
end
%% qu

%%
opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_PixList',PixList,...
    'FPD_RingList',RingList,...
    'FPD_MeanEff',FPD_MeanEff,...
    'FPD_ROIEff',FPD_ROIEff,...
    'FPD_ROIlow',FPD_ROIlow,...
    'FPD_PileUpEff',FPD_PileUpEff};

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type',BKG_Type,...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec,...
    'BKG_RateRingSec',BKG_RateRingSec};

opt_fsd= {...
    'TTFSD',TTFSD,...
    'DTFSD',DTFSD,...
    'HTFSD',HTFSD};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};

opt_corr = {'RadiativeFlag','ON'};

TritiumObject = TBD(...
    opt_corr{:},...
    opt_calc{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_doppler{:},...
    opt_integration{:});
% ---------------------------------------------------------------------- %
end

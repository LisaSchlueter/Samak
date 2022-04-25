function TritiumObject = ref_KSNX_KATRIN_Final2(varargin)



% use MTD from KNM1
MTD = importdata([getenv('SamakPath'),'tritium-data/mat/Knm1/51536.mat']);

% ---------------------------------------------------------------------- %
% WGTS
p = inputParser;
p.addParameter('TDMode','DataTBD',@(x)ischar(x));
p.addParameter('WGTS_CD_MolPerCm2',5*1e17,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2_SubRun','',@(x)isfloat(x));
p.addParameter('WGTS_MolFrac_TT',0.95,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT_SubRun','');
p.addParameter('WGTS_MolFrac_DT',0.011,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_DT_SubRun','');
p.addParameter('WGTS_MolFrac_HT',0.035,@(x)isfloat(x)&& x>0);
p.addParameter('WGTS_MolFrac_HT_SubRun','');
p.addParameter('WGTS_B_T',2.52,@(x)isfloat(x) && x>0);

% Theory
p.addParameter('ISCS','Edep',@(x)ismember(x,{'Aseev','Theory','Edep'}));
p.addParameter('DTFSD','KNM2_0p5eV',@(x)ismember(x,{'OFF','DOSS','BlindingKNM1','SibilleFull','Sibille0p5eV','BlindingKNM2','KNM2_0p5eV'}));
p.addParameter('HTFSD','KNM2_0p5eV',@(x)ismember(x,{'OFF','SAENZ','BlindingKNM1','SibilleFull','Sibille0p5eV','BlindingKNM2','KNM2_0p5eV'}));
p.addParameter('TTFSD','KNM2_0p5eV',@(x)ismember(x,{'OFF','DOSS','SAENZ','BlindingKNM1','SibilleFull','Sibille0p5eV','BlindingKNM2','KNM2_0p5eV'}));
p.addParameter('ELossFlag','KatrinT2A20',@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','KatrinD2','KatrinT2','KatrinT2A20'}));
p.addParameter('DopplerEffectFlag','FSD',@(x)ismember(x,{'OFF','FSD'}));
p.addParameter('AngularTFFlag','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('SynchrotronFlag','ON',@(x)ismember(x,{'OFF','ON'}));

% Binning
p.addParameter('nTeBinningFactor',100,@(x)isfloat(x) && x>0);
p.addParameter('qU',mean(MTD.qU,2),@(x)isfloat(x) && all(x>0));
p.addParameter('qUfrac',mean(MTD.qUfrac,2),@(x)isfloat(x));

% MACE
p.addParameter('MACE_Bmax_T',4.23,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Ba_T',6.3*1e-04,@(x)isfloat(x) && x>0);
p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));
p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('UseParallelRF','ON',@(x)ismember(x,{'OFF','ON'}));

% General
p.addParameter('TimeSec',750*24*60*60,@(x)isfloat(x));
p.addParameter('mnuSq_i',0,@(x)isfloat(x));
p.addParameter('Q_i',18573.7,@(x)isfloat(x));

% FPD
p.addParameter('FPD_MeanEff',0.95,@(x)isfloat(x) && x>0);
p.addParameter('FPD_ROIlow',14,@(x)isfloat(x) && x>0);
p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',0.130);
p.addParameter('BKG_RatePixelSec',''); 
p.addParameter('BKG_RateRingSec','');
p.addParameter('FPD_Segmentation','OFF',@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));
p.addParameter('PixList',1:136);
p.addParameter('RingList',1:13);

p.parse(varargin{:});

WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_CD_MolPerCm2_SubRun = p.Results.WGTS_CD_MolPerCm2_SubRun;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
WGTS_MolFrac_HT_SubRun   = p.Results.WGTS_MolFrac_HT_SubRun;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
WGTS_MolFrac_DT_SubRun   = p.Results.WGTS_MolFrac_DT_SubRun;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;
WGTS_MolFrac_TT_SubRun   = p.Results.WGTS_MolFrac_TT_SubRun;
WGTS_B_T                 = p.Results.WGTS_B_T;

% Parameters: Calculation precision
TDMode                   = p.Results.TDMode;
nTeBinningFactor         = p.Results.nTeBinningFactor;
qU                       = p.Results.qU;
qUfrac                   = p.Results.qUfrac;

%MACE 
MACE_Bmax_T              = p.Results.MACE_Bmax_T;
MACE_Ba_T                = p.Results.MACE_Ba_T;
KTFFlag                  = p.Results.KTFFlag;
recomputeRF              = p.Results.recomputeRF;
UseParallelRF            = p.Results.UseParallelRF;

%Theory
ISCS                     = p.Results.ISCS;
ELossFlag                = p.Results.ELossFlag;
DTFSD                    = p.Results.DTFSD;
HTFSD                    = p.Results.HTFSD;
TTFSD                    = p.Results.TTFSD;
DopplerEffectFlag        = p.Results.DopplerEffectFlag;
AngularTFFlag            = p.Results.AngularTFFlag;
SynchrotronFlag          = p.Results.SynchrotronFlag;

%FPD
FPD_Segmentation         = p.Results.FPD_Segmentation;
FPD_MeanEff              = p.Results.FPD_MeanEff;
FPD_ROIlow               = p.Results.FPD_ROIlow;
FPD_ROIEff               = p.Results.FPD_ROIEff;
FPD_PileUpEff            = p.Results.FPD_PileUpEff;
BKG_Flag                 = p.Results.BKG_Flag;
BKG_Type                 = p.Results.BKG_Type;
BKG_RateAllFPDSec        = p.Results.BKG_RateAllFPDSec;
PixList                  = p.Results.PixList;
RingList                 = p.Results.RingList;

%General
mnuSq_i                  = p.Results.mnuSq_i;
Q_i                      = p.Results.Q_i;
TimeSec                  = p.Results.TimeSec;

if isempty(WGTS_MolFrac_TT_SubRun)
    WGTS_MolFrac_TT_SubRun = repmat(WGTS_MolFrac_TT,size(qU));   
end
if isempty(WGTS_MolFrac_DT_SubRun)
    WGTS_MolFrac_DT_SubRun = repmat(WGTS_MolFrac_DT,size(qU));   
end
if isempty(WGTS_MolFrac_HT_SubRun)
    WGTS_MolFrac_HT_SubRun = repmat(WGTS_MolFrac_HT,size(qU));   
end

if isempty(WGTS_CD_MolPerCm2_SubRun)
    WGTS_CD_MolPerCm2_SubRun = repmat(WGTS_CD_MolPerCm2,size(qU));   
end

% ---------------------------------------------------------------------- %
% Create Tritium spectrum object
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD','DataDriven',...
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
    'KTFFlag',KTFFlag,...
    'AngularTFFlag',AngularTFFlag,...
    'SynchrotronFlag',SynchrotronFlag};

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
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};

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
end



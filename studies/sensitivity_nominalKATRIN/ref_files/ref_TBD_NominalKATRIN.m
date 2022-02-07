function TritiumObject = ref_TBD_NominalKATRIN(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration File for sensitivity studies 
% KATRIN nominal settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parser
p = inputParser;

% TBDDS 
p.addParameter('mnuSq_i',0,@(x)isfloat(x));
p.addParameter('Q_i',18574,@(x)isfloat(x)); % effective Endpoint

% WGTS
p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT',0.95,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_DT',0.04,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_HT',0.01,@(x)isfloat(x)&& x>0);

% Binning
p.addParameter('nTeBinningFactor',100,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('TD','DR30',@(x) ischar(x));
p.addParameter('WGTS_B_T',3.6,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Bmax_T',6,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Ba_T',3e-4,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Ba_Setting','2.7G',@(x) ischar(x)); %Do NOT use Pixelwise Ba corrections
p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));
p.addParameter('ELossFlag','Abdurashitov',@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','CW_G2LT','KatrinD2'}));

% TD Reading Mode
p.addParameter('TDMode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS','WGTS100K','BlindingKNM1'}));
p.addParameter('HTFSD','SAENZ',@(x)ismember(x,{'OFF','SAENZ','WGTS100K','BlindingKNM1'}));
p.addParameter('TTFSD','SAENZ',@(x)ismember(x,{'OFF','DOSS','SAENZ','WGTS100K','BlindingKNM1'}));

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));

% FPD
p.addParameter('c',1:119,@(x)all(isfloat(x)) && all(x)>0); 
p.addParameter('FPD_MeanEff',0.95,@(x)isfloat(x) && x>0); 
p.addParameter('FPD_ROIlow',14,@(x)isfloat(x) && x>0);  % sets the FPD coverage
p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',0.01,@(x)isfloat(x) && x>=0);
p.addParameter('FPD_Segmentation','OFF',@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));

% Transmission Function
p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));
p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('UseParallelRF','ON',@(x)ismember(x,{'OFF','ON'}));

% Multipixel
p.addParameter('FPD_PixList',1,@(x)isfloat(x) && all(x>=0));

% Ring
p.addParameter('FPD_RingList',1,@(x)isfloat(x) && all(x>=0));


p.parse(varargin{:});

%RunSummary Parameter
TimeSec                  = p.Results.TimeSec;
WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;
TD                       = p.Results.TD;

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
TDMode            = p.Results.TDMode;
WGTS_B_T          = p.Results.WGTS_B_T;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
MACE_Ba_T         = p.Results.MACE_Ba_T; 
MACE_Ba_Setting   = p.Results.MACE_Ba_Setting;
ISCS              = p.Results.ISCS;
ELossFlag         = p.Results.ELossFlag;

FPD_PixList       = p.Results.FPD_PixList;
FPD_MeanEff       = p.Results.FPD_MeanEff;
FPD_ROIlow        = p.Results.FPD_ROIlow;
FPD_ROIEff        = p.Results.FPD_ROIEff;
FPD_PileUpEff     = p.Results.FPD_PileUpEff;
BKG_Flag          = p.Results.BKG_Flag;
BKG_Type          = p.Results.BKG_Type;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;
Q_i               = p.Results.Q_i;
FPD_Segmentation  = p.Results.FPD_Segmentation;

% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD           = p.Results.DTFSD;
HTFSD           = p.Results.HTFSD;
TTFSD           = p.Results.TTFSD;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;
recomputeRF     = p.Results.recomputeRF;
UseParallelRF   = p.Results.UseParallelRF;

% Ring
FPD_RingList        = p.Results.FPD_RingList;

%%%%%%%%%%%%%%% END parser



% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i,...
    'Q_i',Q_i};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30,...
    'ISCS',ISCS,...
    'ELossFlag',ELossFlag,...
    'recomputeRF',recomputeRF,...
    'UseParallelRF',UseParallelRF};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T,...
    'MACE_Ba_Setting',MACE_Ba_Setting};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_PixList',FPD_PixList,...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_RingList',FPD_RingList,...
    'FPD_MeanEff',FPD_MeanEff,...
    'FPD_ROIEff',FPD_ROIEff,...
    'FPD_PileUpEff',FPD_PileUpEff,...
    'FPD_ROIlow',FPD_ROIlow};

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

opt_corr = {...
    'RadiativeFlag','ON',...
    'RecoilWmVmAFlag','ON',...
    'FiniteExtChargeFlag','ON',...
    'EEexchangeFlag','ON',...    
    'ScreeningFlag','OFF',...
    'WintFiniteSizeFlag','ON',...
    'RecoilCoulombFlag','ON'};

% Tritium spectrum definition
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

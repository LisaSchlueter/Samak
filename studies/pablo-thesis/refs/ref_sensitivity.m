function TritiumObject = ref_sensitivity(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration File for the revised sensitivity studies for my master
% thesis
% P. I. M. G. 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parser
p = inputParser;

% Parameter from RunSummary
p.addParameter('TimeSec',3*365.25*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT',0.95,@(x)isfloat(x) && x>0);
% p.addParameter('WGTS_MolFrac_TT_SubRun',D.WGTS_MolFrac_TT_SubRun);
p.addParameter('WGTS_MolFrac_DT',0.04,@(x)isfloat(x) && x>0);
% p.addParameter('WGTS_MolFrac_DT_SubRun',D.WGTS_MolFrac_DT_SubRun);
p.addParameter('WGTS_MolFrac_HT',0.01,@(x)isfloat(x)&& x>0);
% p.addParameter('WGTS_MolFrac_HT_SubRun',D.WGTS_MolFrac_HT_SubRun);
p.addParameter('TD','Flat60',@(x) ischar(x));

% Binning
p.addParameter('nTeBinningFactor',50,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('WGTS_B_T',3.6,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Bmax_T',6,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Ba_T',3e-4,@(x)isfloat(x) && x>0);
p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));

% TD Reading Mode
p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));
p.addParameter('mnuSq_i',0,@(x)isfloat(x));

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('HTFSD','SAENZ',@(x)ismember(x,{'OFF','SAENZ'}));
p.addParameter('TTFSD','SAENZ',@(x)ismember(x,{'OFF','DOSS','SAENZ'}));

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));

% FPD
p.addParameter('FPD_MeanEff',0.95,@(x)isfloat(x) && x>0);
p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',0.01,@(x)isfloat(x) && x>=0);
p.addParameter('FPD_Segmentation','OFF',@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));

% Transmission Function
p.addParameter('KTFFlag','Compute',@(x)ismember(x,{'OFF','MACE','Compute'}));
p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('UseParallelRF','ON',@(x)ismember(x,{'OFF','ON'}));

% Multipixel
p.addParameter('FPD_Pixel',1,@(x)isfloat(x) && all(x>=0));

% Ring
p.addParameter('FPD_Ring',1,@(x)isfloat(x) && all(x>=0));


p.parse(varargin{:});

%RunSummary Parameter
TimeSec                  = p.Results.TimeSec;
WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
%WGTS_MolFrac_HT_SubRun   = p.Results.WGTS_MolFrac_HT_SubRun;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
%WGTS_MolFrac_DT_SubRun   = p.Results.WGTS_MolFrac_DT_SubRun;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;
%WGTS_MolFrac_TT_SubRun   = p.Results.WGTS_MolFrac_TT_SubRun;
TD                       = p.Results.TD;

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
Mode              = p.Results.Mode;
WGTS_B_T          = p.Results.WGTS_B_T;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
MACE_Ba_T         = p.Results.MACE_Ba_T; 
ISCS              = p.Results.ISCS;

FPD_MeanEff       = p.Results.FPD_MeanEff;
FPD_ROIEff        = p.Results.FPD_ROIEff;
FPD_PileUpEff     = p.Results.FPD_PileUpEff;
BKG_Flag          = p.Results.BKG_Flag;
BKG_Type          = p.Results.BKG_Type;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;
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

% Multipixel
FPD_Pixel       = p.Results.FPD_Pixel;

% Ring
FPD_Ring        = p.Results.FPD_Ring;

%%%%%%%%%%%%%%% END parser



% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'Mode',Mode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

%    'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...    

%    'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...

%    'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
%    'WGTS_MolFracRelErr_TT',std(WGTS_MolFrac_TT_SubRun)/WGTS_MolFrac_TT,...
%    'WGTS_MolFracRelErr_DT',std(WGTS_MolFrac_DT_SubRun)/WGTS_MolFrac_DT,...
%    'WGTS_MolFracRelErr_HT',std(WGTS_MolFrac_HT_SubRun)/WGTS_MolFrac_HT,...
opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30,...
    'ISCS',ISCS,...
    'recomputeRF',recomputeRF,...
    'UseParallelRF',UseParallelRF};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_Pixel',FPD_Pixel,...
    'FPD_Ring',FPD_Ring,...
    'FPD_MeanEff',FPD_MeanEff,...
    'FPD_ROIEff',FPD_ROIEff,...
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

function myobj = ref_katrin(varargin)
%
% Configurate KATRIN experiment
%
% Th. Lasserre, December 2017
%

% Parser
p = inputParser;

% Te Binning
p.addParameter('nTeBinningFactor',10,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('TDMode', 'Sim', @(x)ismember(x,{'DataKr', 'DataTBD', 'Sim', 'Read'}));
p.addParameter('TDProfile', 'flat', @(x)ismember(x,{'flat'}));
p.addParameter('qUmin',18545,@(x)isfloat(x) && x>0);
p.addParameter('qUmax',18580,@(x)isfloat(x) && x>0);
p.addParameter('qUStepSize',1,@(x)isfloat(x) && x>0);
p.addParameter('TD','Run40257');
       
p.addParameter('TimeSec',3*365*86400,@(x)isfloat(x) && x>0);
p.addParameter('PS_Wein93','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_RateAllFPDSec',300e-3,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_i',0,@(x)isfloat(x));


% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('TTFSD','SAENZ',@(x)ismember(x,{'OFF','DOSS','SAENZ','SAENZNOEE','ROLL'}));
p.addParameter('DTFSD','OFF',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('HTFSD','OFF',@(x)ismember(x,{'OFF','SAENZ'}));

p.addParameter('WGTS_MolFrac_TT',0.95,@(x)isfloat(x));
p.addParameter('WGTS_CD_MolPerCm2',5e17,@(x)isfloat(x));
    
% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{ 'OFF','Voigt','Conv','realConv'}));

% Transmission Function
p.addParameter('KTFFlag','WGTSMACE');

% Flags for Theoretical Corrections 
p.addParameter('RadiativeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('ScreeningFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('RecoilWmVmAFlag','OFF',@(x)ismember(x,{'OFF','ON'}));


p.parse(varargin{:});

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
TDMode            = p.Results.TDMode;
TDProfile         = p.Results.TDProfile;
qUmin             = p.Results.qUmin;
qUmax             = p.Results.qUmax;
qUStepSize        = p.Results.qUStepSize;
       
TimeSec           = p.Results.TimeSec;
TD                = p.Results.TD;
PS_Wein93         = p.Results.PS_Wein93;

BKG_Flag          = p.Results.BKG_Flag;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;


% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
TTFSD           = p.Results.TTFSD;
DTFSD           = p.Results.DTFSD;
HTFSD           = p.Results.HTFSD;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;

WGTS_MolFrac_TT = p.Results.WGTS_MolFrac_TT;
WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;

% Theoretical Corrections
RadiativeFlag   = p.Results.RadiativeFlag;
ScreeningFlag   = p.Results.ScreeningFlag;
RecoilWmVmAFlag = p.Results.RecoilWmVmAFlag;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_B_T',3.6};

opt_mace = {...
    'MACE_Bmax_T',6,...
    'MACE_Ba_T',6e-4};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_Segmentation','OFF','FPD_MeanEff',0.95};

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};

opt_fsd= {...
    'TTFSD',TTFSD,...
    'DTFSD',DTFSD,...
    'HTFSD',HTFSD;
    };

opt_theocorr = {...
    'ScreeningFlag',ScreeningFlag,...
    'FiniteExtChargeFlag','OFF',...
    'EEexchangeFlag','OFF',...
    'RecoilCoulombFlag','OFF',...
    'RadiativeFlag',RadiativeFlag,...
    'RecoilWmVmAFlag',RecoilWmVmAFlag,...
    'WintFiniteSizeFlag','OFF',...
    'PS_Wein93',PS_Wein93,...
    };

opt_integration = {'IStype','SIMPFAST'};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

% Tritium spectrum definition
myobj = TBD(...
    opt_calc{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_theocorr{:},...
    opt_doppler{:},...
    opt_integration{:});
end


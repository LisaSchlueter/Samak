function myobj = InitKATRINDC(varargin)
%

%

% Parser
p = inputParser;

% Binning
p.addParameter('nTeBinningFactor',100,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
p.addParameter('TimeYear',3,@(x)isfloat(x) && x>0);
p.addParameter('TD','DR30',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn',...
    'Flat5','Flat20','Flat30','Flat50','Flat60','Flat100'}));
p.addParameter('PS_Wein93','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_RateAllFPDSec',100e-3,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_i',0,@(x)isfloat(x));
p.addParameter('Mode','Sim',@(x)ismember(x,{'DataTBD', 'Sim'}));


% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('TTFSD','OFF',@(x)ismember(x,{'OFF','DOSS','SAENZ','SAENZNOEE','ROLL'}));
p.addParameter('DTFSD','OFF',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('HTFSD','OFF',@(x)ismember(x,{'OFF','SAENZ'}));
p.addParameter('T2purity',0.955,@(x)isfloat(x) && x>0);

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','ON'}));

% Transmission Function
p.addParameter('KTFFlag','MACE+WGTSIS',@(x)ismember(x,{'OFF','SSC','SSC2',...
    'SSCMW','MACE','MACEMartin','MACE+WGTSIS','TL','SG','Kle','MACE+RIPPLE',...
    'SSCW_DCperqU','SSCW_DCOfficial','SSCW_DC18545','SSC_BC'}));

% Flags for Theoretical Corrections 
p.addParameter('RadiativeFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('ScreeningFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('RecoilWmVmAFlag','OFF',@(x)ismember(x,{'OFF','ON'}));


p.parse(varargin{:});

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
Normalization     = p.Results.Normalization;
TimeYear          = p.Results.TimeYear;
TD                = p.Results.TD;
PS_Wein93         = p.Results.PS_Wein93;
Mode              = p.Results.Mode;


BKG_Flag          = p.Results.BKG_Flag;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;


% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
TTFSD           = p.Results.TTFSD;
DTFSD           = p.Results.DTFSD;
HTFSD           = p.Results.HTFSD;
T2purity        = p.Results.T2purity;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;

% Theoretical Corrections
RadiativeFlag = p.Results.RadiativeFlag;
ScreeningFlag = p.Results.ScreeningFlag;
RecoilWmVmAFlag = p.Results.RecoilWmVmAFlag;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TD',TD,...
    'TimeYear',TimeYear,...
    'mnuSq_i',mnuSq_i,'Mode',Mode};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.632455255110580,...
    'WGTS_Tp',T2purity,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',5e17,...
    'WGTS_B_T',3.6};

opt_mace = {...
    'MACE_Bmax_T',6,...
    'MACE_Ba_T',3e-4,...
    'MACE_R_eV',0.93};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_Segmentation','OFF'};

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
    opt_doppler{:},'recomputeKernel','OFF');
end


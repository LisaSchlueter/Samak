function myobj = InitKATRINE0_allPixels(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION FILE
% Endpoint sensitivity study 
% 148 pixels
% Configuration for First Tritium May

% P. Morales 2018
% Last update 16/03/2018


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parser
p = inputParser;

% Binning
p.addParameter('nTeBinningFactor',10,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('TimeSec',2*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('TD','Flat30',@(x)ismember(x,{'DR20','DR30','DR40','DR50',...
    'Flat20','Flat30','Flat50','Flat60','Flat100'}));
p.addParameter('mnuSq_i',0,@(x)isfloat(x));
p.addParameter('qUmin',18545,@(x)isfloat(x) && x>0.01);
p.addParameter('qUmax',18580,@(x)isfloat(x) && x>0.01);
p.addParameter('qUStepSize',1,@(x)isfloat(x) && x>0);
p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));

% Background
p.addParameter('BKG_Flag','XmasData',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_RateAllFPDSec',0e-3,@(x)isfloat(x) && x>=0)

% number of pixels
p.addParameter('nPixels',148,@(x)isfloat(x) && x>0);

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('T2purity',0.01,@(x)isfloat(x) && x>0);

% Flag for Doppler Effect
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','ON'}));

% Transmission Function
p.addParameter('KTFFlag','MACE+WGTSIS',@(x)ismember(x,{'OFF','SSC',...
    'SSCMW','MACE','MACEMartin','MACE+WGTSIS','TL','SG','Kle','MACE+RIPPLE',...
    'SSCW_DCperqU','SSCW_DCOfficial','SSCW_DC18545','SSC_BC','Compute'}));

p.parse(varargin{:});

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
TD                = p.Results.TD;
Mode              = p.Results.Mode;
qUmin             = p.Results.qUmin;
qUmax             = p.Results.qUmax;
qUStepSize        = p.Results.qUStepSize;
BKG_Flag          = p.Results.BKG_Flag;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;
TimeSec           = p.Results.TimeSec;

% Background / FPD

% Pixels
nPixels         = p.Results.nPixels;

% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD           = p.Results.DTFSD;
T2purity        = p.Results.T2purity;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'Mode',Mode,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

switch Mode
    case 'Read'
        opt_TD = {...
            'TD',TD};
    case 'Sim'
        opt_TD = {...
            'qUmin',qUmin,...
            'qUmax',qUmax,...
            'qUStepSize',qUStepSize};
end

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.630570414481244,...
    'WGTS_Tp',T2purity,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',2.5e17,...
    'WGTS_B_T',2.52,...
    'WGTS_Temp',30};

opt_mace = {...
    'MACE_Bmax_T',4.2,...
    'MACE_Ba_T',6e-4};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_Segmentation','PIXEL',...
    'FPD_Pixel',0,...
    'nPixels',nPixels};

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};

opt_fsd= {...
    'TTFSD','OFF',...
    'DTFSD',DTFSD,...
    'HTFSD','OFF'};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};

% Tritium spectrum definition
myobj = TBD(...
    opt_calc{:},...
    opt_katrin{:},...
    opt_TD{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_doppler{:},...
    opt_integration{:});

   %'WGTS_MolFrac_DT',1
end


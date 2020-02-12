function myobj = krl332minuit_init(varargin)
%
% Configurate KATRIN experiment
% Krypton Mode
%
% Th. Lasserre, August 2017
%

% Parser
p = inputParser;

% Miscellaneous
p.addParameter('Display','OFF',@(x)ismember(x,{'ON','OFF'}));

% Calculation precision

p.addParameter('nTeBinningFactor',30,@(x)isfloat(x) && x>0);
p.addParameter('Tolerance', 1e-2,@(x)isfloat(x) && x>0);
         
% Krypton Mode
p.addParameter('L3_32_Phi0_i',1e1,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('nPixels',148,@(x)isfloat(x) && x>0);
p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TimeSec',1,@(x)isfloat(x) && x>0);
p.addParameter('TD','MoSFitterJuly2017',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}));
p.addParameter('BKG_RateAllFPDSec',1e-3,@(x)isfloat(x) && x>0);
% FPD
p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','PIXEL'}));
% MAC E
p.addParameter('HVRipples','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
% Doppler
p.addParameter('DopplerEffectFlag','Voigt',@(x)ismember(x,{'OFF','matConv','Voigt','numConv'}));
p.addParameter('WGTS_Temp',100,@(x)isfloat(x) && x>0);
p.addParameter('DE_sigma',0.5,@(x)isfloat(x) && x>0); 

p.parse(varargin{:});

% Miscellaneous
Display = p.Results.Display;

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;
Tolerance        = p.Results.Tolerance;

% Krypton Mode
L3_32_Phi0_i      = p.Results.L3_32_Phi0_i;

% Parameters: KATRIN GENERAL SETTINGS
nPixels           = p.Results.nPixels;
Normalization     = p.Results.Normalization;
CPS               = p.Results.CPS;
TimeSec          = p.Results.TimeSec;
TD                = p.Results.TD;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
% FPD
FPD_Segmentation  = p.Results.FPD_Segmentation;
% MACE
HVRipples         = p.Results.HVRipples;
HVRipplesP2PV     = p.Results.HVRipplesP2PV;
% Doppler
DopplerEffectFlag = p.Results.DopplerEffectFlag;
WGTS_Temp         = p.Results.WGTS_Temp;
DE_sigma          = p.Results.DE_sigma;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor,...
    'Tolerance',Tolerance,...
    'CPS',CPS,...
    'nPixels',nPixels};

opt_krypton = {...
    'L3_32_Phi0_i',L3_32_Phi0_i,...
    };

opt_katrin = {...
    'TD',TD,...
    'Normalization',Normalization,...
    'TimeSec',TimeSec};

opt_wgts = {...
            'WGTS_CosMaxAAngle',0.6324,...
            'WGTS_Tp',0.95,...
            'WGTS_DTHTr',0.1,...
            'WGTS_FTR_cm',4.1,...
            'WGTS_CD_MolPerCm2',5e17,...
            'WGTS_B_T',3.6};

switch TD
    case {'MoSFitterJuly2017'}
        Eline = 30472.300;
end

opt_mace = {...
    'MACE_Bmax_T',6,...
    'MACE_Ba_T',3e-04,...
    'MACE_R_eV',3e-04/6*Eline,...
    'HVRipples','OFF'};

opt_wgtsmace = {...
    'KTFFlag','MACEMartin'};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};


opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag,...
            'recomputeKernel','OFF','DE_sigma',DE_sigma,'WGTS_Temp',WGTS_Temp,'Eplus',4,'Eminus',4};


% Tritium spectrum definition
myobj = Kr(...
    opt_krypton{:},...
    opt_calc{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_bkg{:},...
    opt_doppler{:});

% myobj.ComputeKrDS;  myobj.PlotKrDS('fig',1);
% myobj.ComputeKrIS;  myobj.PlotKrIS('fig',2);

switch Display
    case 'ON'
        myobj.DisplayKrInfo
end

end


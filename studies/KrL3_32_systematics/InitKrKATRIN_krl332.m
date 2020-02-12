function myobj = InitKrKATRIN_krl332(varargin)
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
p.addParameter('Mode', 'DataKr');
% Calculation precision
p.addParameter('nTeBinningFactor',40,@(x)isfloat(x) && x>0);
p.addParameter('Tolerance',1e-2,@(x)isfloat(x) && x>0);       
% Krypton Mode
p.addParameter('L3_32_Phi0_i',1e1,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('nPixels',148,@(x)isfloat(x) && x>0);
p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
p.addParameter('CPS','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','KrL3_32',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}));
p.addParameter('BKG_RateAllFPDSec',1e-3,@(x)isfloat(x) && x>0);
p.addParameter('BKG_RateSecallPixels',5,@(x)isfloat(x) && x>0);
% FPD
p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'OFF','RING','PIXEL'}));
% MAC E
p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.44,@(x)isfloat(x) && x>0); % V or eV
% Doppler
p.addParameter('DopplerEffectFlag','Voigt',@(x)ismember(x,{'OFF','matConv','Voigt','numConv'}));
p.addParameter('DE_sigma',0.058,@(x)isfloat(x) && x>0); 

p.parse(varargin{:});

% Miscellaneous
Display = p.Results.Display;
Mode    = p.Results.Mode;
% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;
Tolerance        = p.Results.Tolerance;
% Krypton Mode
L3_32_Phi0_i      = p.Results.L3_32_Phi0_i;

% Parameters: KATRIN GENERAL SETTINGS
nPixels           = p.Results.nPixels;
Normalization     = p.Results.Normalization;
CPS               = p.Results.CPS;
TD                = p.Results.TD;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
BKG_RateSecallPixels = p.Results.BKG_RateSecallPixels;
% FPD
FPD_Segmentation  = p.Results.FPD_Segmentation;
% MACE
HVRipples         = p.Results.HVRipples;
HVRipplesP2PV     = p.Results.HVRipplesP2PV;
% Doppler
DopplerEffectFlag = p.Results.DopplerEffectFlag;
DE_sigma          = p.Results.DE_sigma;



% Default
opt_calc = {...
    'Mode', Mode,...
    'nTeBinningFactor',nTeBinningFactor,...
    'Tolerance',Tolerance,...
    'CPS',CPS,...
    'nPixels',nPixels};

opt_krypton = {...
    'L3_32_Phi0_i',L3_32_Phi0_i,...
    };

opt_katrin = {...
    'TD',TD,...
    'Normalization',Normalization};

opt_wgts = {...
            'WGTS_CosMaxAAngle',0.8901,...%acceptance angle (radian)
            'WGTS_Tp',0.95,...
            'WGTS_DTHTr',0.1,...
            'WGTS_FTR_cm',4.1,...
            'WGTS_CD_MolPerCm2',5e17,...&column density
            'WGTS_B_T',2.5}; %B-field source

switch TD
    case {'MoSFitterJuly2017'}
        Eline = 30473;
    case{'KrL3_32'}
        Eline  = 30472.5;
end

opt_mace = {...
    'MACE_Bmax_T',4.2,...%maximal B-field
    'MACE_Ba_T',2.7e-04,...%minimal B-field
    'MACE_R_eV',2.7e-04/4.2*Eline,...&Energy Resolution
    'HVRipples', HVRipples, ...
    'HVRipplesP2PV', HVRipplesP2PV};

opt_wgtsmace = {...
    'KTFFlag','MACEMartin'};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec,...
    'BKG_RateSecallPixels', BKG_RateSecallPixels};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag, 'DE_sigma', DE_sigma};

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


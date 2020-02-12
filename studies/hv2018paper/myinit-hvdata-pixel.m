function myobj = init-hvdata-pixel(varargin)
%
% Configurate KATRIN experiment
% Krypton L3-32 HV Paper Data
%
% Th. Lasserre
% Last Update : Jan 2018
%

% Parser
p = inputParser;

% Miscellaneous
p.addParameter('Display','OFF',@(x)ismember(x,{'ON','OFF'}));

% Calculation precision
p.addParameter('nTeBinningFactor',50,@(x)isfloat(x) && x>0);
         
%
p.addParameter('nPixels',40,@(x)isfloat(x) && x>0);

% Krypton Mode
p.addParameter('L3_32_Phi0_i',60,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TimeYear',1/365.25/86400*30*43,@(x)isfloat(x) && x>0);
p.addParameter('TD','DummyKrdc1',@(x)ismember(x,{'DummyKrdc1','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}));
p.addParameter('BKG_RateAllFPDSec', 10,@(x)isfloat(x) && x>0);
p.addParameter('BKG_RateRingSec',   10,@(x)isfloat(x) && x>0);
p.addParameter('BKG_RatePixelSec',  10,@(x)isfloat(x) && x>0);
% FPD
p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('FPD_Pixel',1,@(x)isfloat(x) && x>-1); % 0=All Pixels
p.addParameter('FPD_Ring',1,@(x)isfloat(x) && x>-1);  % 0=All Rings
% MAC E
p.addParameter('HVRipples','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
p.addParameter('KTFFlag', 'MACER',@(x)ismember(x,{'OFF','SSC','MACE','MACER','MACE+WGTSIS','TL','SG','Kle'}));
% Doppler
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','Conv','Voigt','realConv'}));
p.addParameter('WGTS_Temp',100,@(x)isfloat(x) && x>0);
p.addParameter('DE_sigma',0.058,@(x)isfloat(x) && x>0); % 0.134429846285963

p.parse(varargin{:});

% Miscellaneous
Display = p.Results.Display;

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

%
nPixels            = p.Results.nPixels;

% Krypton Mode
L3_32_Phi0_i      = p.Results.L3_32_Phi0_i;

% Parameters: KATRIN GENERAL SETTINGS
Normalization     = p.Results.Normalization;
CPS               = p.Results.CPS;
TimeYear          = p.Results.TimeYear;
TD                = p.Results.TD;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
BKG_RateRingSec   = p.Results.BKG_RateRingSec;
BKG_RatePixelSec  = p.Results.BKG_RatePixelSec;
% FPD
FPD_Segmentation  = p.Results.FPD_Segmentation;
FPD_Pixel         = p.Results.FPD_Pixel;
FPD_Ring          = p.Results.FPD_Ring;
% MACE
HVRipples         = p.Results.HVRipples;
HVRipplesP2PV     = p.Results.HVRipplesP2PV;
KTFFlag           = p.Results.KTFFlag;
% Doppler
DopplerEffectFlag = p.Results.DopplerEffectFlag;
DE_sigma          = p.Results.DE_sigma;


% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor,...
    'CPS',CPS,...
    'nPixels',nPixels};

opt_krypton = {...
    'L3_32_Phi0_i',L3_32_Phi0_i};

opt_katrin = {...
    'TD',TD,...
    'Normalization',Normalization,...
    'TimeYear',TimeYear};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.6324,...
    'WGTS_Tp',0.95,...
    'WGTS_DTHTr',0.1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',5e17,...
    'WGTS_B_T',1.8};

Eline = 30472.5;

opt_mace = {...
    'MACE_Bmax_T',4.199,...
    'MACE_Ba_T',2.71e-04,...
    'MACE_R_eV',2.7e-04/4.199*Eline,...
    'HVRipples',HVRipples,...
    'HVRipplesP2PV',HVRipplesP2PV};

opt_wgtsmace = {...
    'KTFFlag', KTFFlag};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_Pixel',FPD_Pixel};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec,...
    'BKG_RateRingSec',BKG_RateRingSec,...
    'BKG_RatePixelSec',BKG_RatePixelSec};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag,...
    'recomputeKernel','OFF','DE_sigma',0.058,'Eplus',4,'Eminus',4};

end


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


function myobj = InitKrKATRIN_AllPixels(varargin)
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
p.addParameter('nTeBinningFactor',19,@(x)isfloat(x) && x>0);
         
% 
p.addParameter('nPixels',148,@(x)isfloat(x) && x>0);

% Krypton Mode
p.addParameter('L3_32_Phi0_i',1e1,@(x)isfloat(x) && x>0);
p.addParameter('K32_Phi0_i',1e1,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
p.addParameter('CPS','ON',@(x)ismember(x,{'ON','OFF'}));
%p.addParameter('TimeYear',1/10000,@(x)isfloat(x) && x>0);
p.addParameter('TimeYear',1/1000,@(x)isfloat(x) && x>0);
p.addParameter('TD','KrL3_32',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}));
p.addParameter('BKG_RateAllFPDSec', 1000e-3,@(x)isfloat(x) && x>0);
p.addParameter('BKG_RateRingSec',1000e-3,@(x)isfloat(x) && x>0);
p.addParameter('BKG_RatePixelSec',1000e-3,@(x)isfloat(x) && x>0);
% FPD
p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('FPD_Pixel',1,@(x)isfloat(x) && x>-1); % 0=All Pixels
p.addParameter('FPD_Ring',1,@(x)isfloat(x) && x>-1);  % 0=All Rings
% MAC E
p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
% Doppler Effect
p.addParameter('DopplerEffectFlag','Conv',@(x)ismember(x,{'OFF','Conv','Voigt','realConv'}));
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
K32_Phi0_i        = p.Results.K32_Phi0_i;

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
% Doppler
DopplerEffectFlag = p.Results.DopplerEffectFlag;
WGTS_Temp         = p.Results.WGTS_Temp;
DE_sigma          = p.Results.DE_sigma;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor,...
    'CPS',CPS,...
    'nPixels',nPixels};

opt_krypton = {...
    'L3_32_Phi0_i',L3_32_Phi0_i,...
    'K32_Phi0_i',K32_Phi0_i,...
    };

opt_katrin = {...
    'TD',TD,...
    'Normalization',Normalization,...
    'TimeYear',TimeYear};

switch TD
    case {'KrL3_32','KrL3_32_HS','KrK32'}
        opt_wgts = {...
            'WGTS_CosMaxAAngle',0.6324,...
            'WGTS_Tp',0.95,...
            'WGTS_DTHTr',0.1,...
            'WGTS_FTR_cm',4.1,...
            'WGTS_CD_MolPerCm2',5e17,...
            'WGTS_B_T',2.52};
    case 'MoSFitterJuly2017'
        opt_wgts = {...
            'WGTS_CosMaxAAngle',0.6324,...
            'WGTS_Tp',0.95,...
            'WGTS_DTHTr',0.1,...
            'WGTS_FTR_cm',4.1,...
            'WGTS_CD_MolPerCm2',5e17,...
            'WGTS_B_T',3.6};
end

switch TD
    case 'KrK32'
        Eline = 17825;
    case {'KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}
        Eline = 30473;
end

switch TD
    case {'KrL3_32','KrL3_32_HS','KrK32'}
        opt_mace = {...
            'MACE_Bmax_T',4.199,...
            'MACE_Ba_T',2.7e-04,...
            'MACE_R_eV',2.7e-04/4.199*Eline,...
            'HVRipples','ON',...
            'HVRipplesP2PV',0.52};
    case 'MoSFitterJuly2017'
        opt_mace = {...
            'MACE_Bmax_T',6,...
            'MACE_Ba_T',3e-04,...
            'MACE_R_eV',3e-04/6*Eline,...
            'HVRipples','OFF'};
end

opt_wgtsmace = {...
    'KTFFlag','MACE'};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_Pixel',FPD_Pixel};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec,...
    'BKG_RateRingSec',BKG_RateRingSec,...
    'BKG_RatePixelSec',BKG_RatePixelSec};

switch TD
    case {'KrL3_32','KrL3_32_HS','KrK32'}
        opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag,'DE_sigma',0.058};
    case 'MoSFitterJuly2017'
        opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag,'DE_sigma',0.5};
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


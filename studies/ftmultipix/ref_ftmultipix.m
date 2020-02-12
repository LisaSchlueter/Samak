function TritiumObject = ref_ftmultipix(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION FILE
% Endpoint sensitivity study 
% Multipixe Detector
% Configuration for First Tritium May

% T. Lasserre 2018
% P. I. Morales Guzm�n 2018
% Last update 23/04/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parser
p = inputParser;

% Binning
p.addParameter('nTeBinningFactor',5,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('TimeSec',24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2',0.9*5e17);
p.addParameter('WGTS_B_T',2.52);
p.addParameter('MACE_Bmax_T',4.2);
p.addParameter('MACE_Ba_T',6e-4);

p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));
p.addParameter('TD','FT-TL4');
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',300e-3,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_i',0,@(x)isfloat(x));

p.addParameter('nPixels',1,@(x)isfloat(x) && x>0);

p.addParameter('FPD_Pixel',1,@(x)isfloat(x) && sum(x)>-1)
p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','MULTIPIXEL','SINGLEPIXEL'}));

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS'}));

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));

% Transmission Function
p.addParameter('KTFFlag','Compute',@(x)ismember(x,{'OFF','SSC',...
    'MACE','MACE+WGTSIS','TL','SG','Kle',...
    'SSCW_DCperqU','SSC_BC','Compute'}));

p.addParameter('UseParallelRF','OFF',@(x)ismember(x,{'OFF','ON'}));


p.parse(varargin{:});

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
Mode              = p.Results.Mode;
TD                = p.Results.TD;
WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
WGTS_B_T          = p.Results.WGTS_B_T;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
MACE_Ba_T         = p.Results.MACE_Ba_T; 

BKG_Flag          = p.Results.BKG_Flag;
BKG_Type          = p.Results.BKG_Type;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;
TimeSec           = p.Results.TimeSec;

% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD           = p.Results.DTFSD;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;

% Multipixel
UseParallelRF   = p.Results.UseParallelRF;
FPD_Pixel       = p.Results.FPD_Pixel;
FPD_Segmentation= p.Results.FPD_Segmentation;



% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'Mode',Mode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.630570414481244,...
    'WGTS_MolFrac_TT',0,...
    'WGTS_MolFrac_DT',1.2e-2,...
    'WGTS_MolFrac_HT',0,...
    'WGTS_MolFracRelErr_TT',0,...
    'WGTS_MolFracRelErr_DT',0,...
    'WGTS_MolFracRelErr_HT',0,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_Pixel',FPD_Pixel};

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type',BKG_Type,...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};

opt_fsd= {...
    'TTFSD','OFF',...
    'DTFSD',DTFSD,...
    'HTFSD','OFF'};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};

opt_multipixel = {...
    'UseParallelRF',UseParallelRF};

% Tritium spectrum definition
TritiumObject = TBD(...
    opt_calc{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_doppler{:},...
    opt_integration{:},...
    opt_multipixel{:},...
    'recomputeRF','ON');
end


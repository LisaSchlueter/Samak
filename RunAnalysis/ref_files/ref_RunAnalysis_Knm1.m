
function TritiumObject = ref_RunAnalysis_Knm1(D,varargin)%,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION FILE
% retrieves all information from data file
% only some exceptions: WGTS_B_T, MACE_Bmax_T, FPD_ROIlow, ...
%
% T. Lasserre
% L. Schlueter
% P. I. Morales Guzman
%
% Last update 04/06/2019 (Lisa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

%General KATRIN 
p.addParameter('TimeSec',D.TimeSec,@(x)isfloat(x));
p.addParameter('mnuSq_i',0,@(x)isfloat(x));
p.addParameter('Q_i',18573.7,@(x)isfloat(x)); 

% Binning
p.addParameter('nTeBinningFactor',100,@(x)isfloat(x) && x>0);
p.addParameter('qU',D.qU,@(x)isfloat(x) && all(x>0));
p.addParameter('qUfrac',D.qUfrac,@(x)isfloat(x));
p.addParameter('TDMode','DataTBD',@(x)ischar(x));
p.addParameter('TD','DataDriven',@(x)ischar(x));  %TD-name not used anymore: here just placeolder -> could also be removed

% WGTS
p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2,@(x)isfloat(x) && x>0); 
p.addParameter('WGTS_CD_MolPerCm2_SubRun',D.WGTS_CD_MolPerCm2_SubRun);
p.addParameter('WGTS_MolFrac_TT',D.WGTS_MolFrac_TT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT_SubRun',D.WGTS_MolFrac_TT_SubRun);
p.addParameter('WGTS_MolFrac_DT',D.WGTS_MolFrac_DT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_DT_SubRun',D.WGTS_MolFrac_DT_SubRun);
p.addParameter('WGTS_MolFrac_HT',D.WGTS_MolFrac_HT,@(x)isfloat(x)&& x>0);
p.addParameter('WGTS_MolFrac_HT_SubRun',D.WGTS_MolFrac_HT_SubRun);
p.addParameter('WGTS_B_T',2.52,@(x)isfloat(x) && x>0);

%MACE
p.addParameter('MACE_Bmax_T',4.23,@(x)isfloat(x) && x>0);
p.addParameter('MACE_Ba_T',D.MACE_Ba_T,@(x)isfloat(x) && x>0);
p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));                      % WGTSMACE== normal response function
p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));                                    % if ON: complete RF calculated from scratch
p.addParameter('UseParallelRF','ON',@(x)ismember(x,{'OFF','ON'}));                                   % calculate response function in parallel

% Theory
p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));                                  % inelastic scattering cross section
p.addParameter('DTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','BlindingKNM1'}));                % final state distributions
p.addParameter('HTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','SAENZ','BlindingKNM1'})); 
p.addParameter('TTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','SAENZ','BlindingKNM1'}));
p.addParameter('ELossFlag','KatrinD2',@(x)ismember(x,{'Aseev','Abdurashitov','CW_GLT','KatrinD2'})); % description of energy loss in scattering
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));

% FPD
p.addParameter('FPD_MeanEff',0.9,@(x)isfloat(x) && x>0);%%
p.addParameter('FPD_ROIlow',14,@(x)isfloat(x) && x>0);
p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec','');
p.addParameter('BKG_RatePixelSec','');
p.addParameter('BKG_RateRingSec','');
p.addParameter('FPD_Segmentation','OFF',@(x) ismember(x,{'OFF','SINGLEPIXEL','MULTIPIXEL','RING'}));
p.addParameter('PixList',1:148,@(x)isfloat(x));
p.addParameter('RingList',1:12,@(x)isfloat(x));

p.parse(varargin{:});

%KATRIN General
mnuSq_i                  = p.Results.mnuSq_i;
Q_i                      = p.Results.Q_i;
TimeSec                  = p.Results.TimeSec;

% WGTS
WGTS_B_T                 = p.Results.WGTS_B_T;
WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_CD_MolPerCm2_SubRun = p.Results.WGTS_CD_MolPerCm2_SubRun;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
WGTS_MolFrac_HT_SubRun   = p.Results.WGTS_MolFrac_HT_SubRun;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
WGTS_MolFrac_DT_SubRun   = p.Results.WGTS_MolFrac_DT_SubRun;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;
WGTS_MolFrac_TT_SubRun   = p.Results.WGTS_MolFrac_TT_SubRun;

% Binning
TD                       = p.Results.TD;
TDMode                   = p.Results.TDMode;
qU                       = p.Results.qU;
qUfrac                   = p.Results.qUfrac;
nTeBinningFactor         = p.Results.nTeBinningFactor;   % Parameters: Calculation precision

%FPD
FPD_Segmentation         = p.Results.FPD_Segmentation;
FPD_MeanEff              = p.Results.FPD_MeanEff;
FPD_ROIlow               = p.Results.FPD_ROIlow;
FPD_ROIEff               = p.Results.FPD_ROIEff;
FPD_PileUpEff            = p.Results.FPD_PileUpEff;
BKG_Flag                 = p.Results.BKG_Flag;
BKG_Type                 = p.Results.BKG_Type;
BKG_RateAllFPDSec        = p.Results.BKG_RateAllFPDSec;
BKG_RateRingSec          = p.Results.BKG_RateRingSec;
BKG_RatePixelSec         = p.Results.BKG_RatePixelSec;
PixList                  = p.Results.PixList;
RingList                 = p.Results.RingList;

% Theory
DopplerEffectFlag        = p.Results.DopplerEffectFlag;
DTFSD                    = p.Results.DTFSD;                 % TBD: Flag FSD's
HTFSD                    = p.Results.HTFSD;
TTFSD                    = p.Results.TTFSD;
ISCS                     = p.Results.ISCS;
ELossFlag                = p.Results.ELossFlag;

% MACE
KTFFlag                  = p.Results.KTFFlag;
recomputeRF              = p.Results.recomputeRF;
UseParallelRF            = p.Results.UseParallelRF;
MACE_Bmax_T              = p.Results.MACE_Bmax_T;
MACE_Ba_T                = p.Results.MACE_Ba_T; 

% ---------------------------------------------------------------------- %
% Create Tritium spectrum object 
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD',TD...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i,...
    'Q_i',Q_i,...
    'qU',qU,...
    'qUfrac',qUfrac};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...    
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
    'WGTS_MolFracRelErr_TT',std(WGTS_MolFrac_TT_SubRun)/WGTS_MolFrac_TT,...
    'WGTS_MolFracRelErr_DT',std(WGTS_MolFrac_DT_SubRun)/WGTS_MolFrac_DT,...
    'WGTS_MolFracRelErr_HT',std(WGTS_MolFrac_HT_SubRun)/WGTS_MolFrac_HT,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.5,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_CD_MolPerCm2_SubRun',WGTS_CD_MolPerCm2_SubRun,...    
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30,...
    'ISCS',ISCS,...
    'ELossFlag',ELossFlag,...
    'recomputeRF',recomputeRF,...
    'UseParallelRF',UseParallelRF,...
    'NIS',7}; %Thierry WARNING 13/2/2019 

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T,...
    'MACE_Ba_Setting','Data'}; %do not use pixelmap correction from txt file, use values from runsummary

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};


%% Get Background: Whole FPD (excluding some pixels), per Ring, per Pixel
switch RingMerge
    case 'Default'
        [PixList,RingPixList] = Ring2PixelDefCombi(RingList,PixList);

    case 'None'
        [PixList,RingPixList] = Ring2Pixel(RingList,PixList);
        
    case 'Full'
        [PixList,RingPixList] = Ring2PixelCombi(RingList,PixList);

end

if isempty(BKG_RateAllFPDSec)
    BKG_RateAllFPDSec = numel(PixList)/148*GetBackground('MACE_Ba_T',mean(MACE_Ba_T),...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
end
if isempty(BKG_RateRingSec) && strcmp(FPD_Segmentation,'RING')
    nPixperRing = cellfun(@(x) numel(x),RingPixList); nPixperRing(~ismember(1:13,RingList)) = [];
    BKG_RateRingSec = (nPixperRing./148)'.*GetBackground('MACE_Ba_T',MACE_Ba_T,...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
    BKG_RateRingSec(isnan(BKG_RateRingSec))=0; BKG_RateRingSec(isinf(BKG_RateRingSec))=0;
elseif isempty(BKG_RateRingSec) && ~strcmp(FPD_Segmentation,'RING')
    BKG_RateRingSec = BKG_RateAllFPDSec./numel(RingList).*ones(numel(RingList),1);
end

if isempty(BKG_RatePixelSec) && ismember(FPD_Segmentation,{'SINGLEPIXEL','MULTIPIXEL'})
    BKG_RatePixelSec = 1/148.*GetBackground('MACE_Ba_T',MACE_Ba_T,...
        'WGTS_B_T',WGTS_B_T,...
        'FPD_ROIlow',FPD_ROIlow);
    BKG_RatePixelSec(isnan(BKG_RatePixelSec))=0;  BKG_RatePixelSec(isinf(BKG_RatePixelSec))=0;
end
%%
opt_fpd = {...
    'FPD_Segmentation',FPD_Segmentation,...
    'FPD_PixList',PixList,...
    'FPD_RingList',RingList,...
    'FPD_RingMerge',RingMerge,...
    'FPD_MeanEff',FPD_MeanEff,...
    'FPD_ROIEff',FPD_ROIEff,...
    'FPD_ROIlow',FPD_ROIlow,...
    'FPD_PileUpEff',FPD_PileUpEff};
    

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type',BKG_Type,...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec,...
    'BKG_RateRingSec',BKG_RateRingSec};

opt_fsd= {...
    'TTFSD',TTFSD,...
    'DTFSD',DTFSD,...
    'HTFSD',HTFSD};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};

opt_corr = {'RadiativeFlag','ON'};

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

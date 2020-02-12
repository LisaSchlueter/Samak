function myobj = ref_runsummaries(RunNr,ex2,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION FILE
% scan down at 0.9% DT - column density stable
%
% T. Lasserre 2018
% L.Schlueter
% Last update 25/05/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load RunSummary
D = importdata([num2str(RunNr),ex2,'.mat']); 

% Parser
p = inputParser;
%p.addParameter('RunNr',40263,@(x)isfloat(x) && x>0);
%p.parse(varargin{:});
%RunNr           = p.Results.RunNr;

% Parameter from RunSummary
p.addParameter('TimeSec',D.TimeSec,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2',D.WGTS_CD_MolPerCm2,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2_SubRun',D.WGTS_CD_MolPerCm2_SubRun);
p.addParameter('WGTS_MolFrac_TT', D.WGTS_MolFrac_TT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_TT_SubRun',D.WGTS_MolFrac_TT_SubRun);
p.addParameter('WGTS_MolFrac_DT', D.WGTS_MolFrac_DT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_DT_SubRun',D.WGTS_MolFrac_DT_SubRun);
p.addParameter('WGTS_MolFrac_HT', D.WGTS_MolFrac_HT,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_HT_SubRun',D.WGTS_MolFrac_HT_SubRun);
p.addParameter('TD',['Run',num2str(RunNr),ex2]);

% Binning
p.addParameter('nTeBinningFactor',50,@(x)isfloat(x) && x>0);
% KATRIN GENERAL SETTINGS
p.addParameter('WGTS_B_T',2.52);
p.addParameter('MACE_Bmax_T',4.2);
p.addParameter('MACE_Ba_T',6e-4);
p.addParameter('ISCS','Theory',@(x)ismember(x,{'Aseev','Theory'}));
p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',0.2,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_i',0,@(x)isfloat(x));

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('HTFSD','SAENZ',@(x)ismember(x,{'OFF','SAENZ'}));
p.addParameter('TTFSD','SAENZ',@(x)ismember(x,{'OFF','SAENZ'}));
p.addParameter('recomputeRF','OFF',@(x)ismember(x,{'ON','OFF'}));
% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));


% Transmission Function
p.addParameter('KTFFlag','Compute',@(x)ismember(x,{'OFF','SSC',...
    'MACE','MACE+WGTSIS','TL','SG','Kle',...
    'SSCW_DCperqU','SSC_BC','Compute'}));

p.parse(varargin{:});

%RunSummary Parameter
TimeSec                  = p.Results.TimeSec;
WGTS_CD_MolPerCm2        = p.Results.WGTS_CD_MolPerCm2;
WGTS_CD_MolPerCm2_SubRun = p.Results.WGTS_CD_MolPerCm2_SubRun;
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT; 
WGTS_MolFrac_HT_SubRun   = p.Results.WGTS_MolFrac_HT_SubRun;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT; 
WGTS_MolFrac_DT_SubRun   = p.Results.WGTS_MolFrac_DT_SubRun;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT; 
WGTS_MolFrac_TT_SubRun   = p.Results.WGTS_MolFrac_TT_SubRun;
TD                       = p.Results.TD;

% Parameters: Calculation precision
nTeBinningFactor         = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
Mode                     = p.Results.Mode;
WGTS_B_T                 = p.Results.WGTS_B_T;
MACE_Bmax_T              = p.Results.MACE_Bmax_T;
MACE_Ba_T                = p.Results.MACE_Ba_T; 
ISCS                     = p.Results.ISCS;

BKG_Flag                 = p.Results.BKG_Flag;
BKG_Type                 = p.Results.BKG_Type;
BKG_RateAllFPDSec        = p.Results.BKG_RateAllFPDSec;
mnuSq_i                  = p.Results.mnuSq_i;

% Doppler Effect
DopplerEffectFlag        = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD                    = p.Results.DTFSD;
HTFSD                    = p.Results.HTFSD;
TTFSD                    = p.Results.TTFSD;

% Transmission Function 
KTFFlag                  = p.Results.KTFFlag;
recomputeRF              = p.Results.recomputeRF;
%%%%%%%%%%%%%%% END parser

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'Mode',Mode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_TT_SubRun',WGTS_MolFrac_TT_SubRun,...    
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_DT_SubRun',WGTS_MolFrac_DT_SubRun,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT,...
    'WGTS_MolFrac_HT_SubRun',WGTS_MolFrac_HT_SubRun,...
    'WGTS_MolFracRelErr_TT',std(WGTS_MolFrac_TT_SubRun),...
    'WGTS_MolFracRelErr_DT',std(WGTS_MolFrac_DT_SubRun),...
    'WGTS_MolFracRelErr_HT',std(WGTS_MolFrac_HT_SubRun),...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_CD_MolPerCm2_SubRun',WGTS_CD_MolPerCm2_SubRun,...
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30,...
    'ISCS',ISCS,...
    'recomputeRF',recomputeRF};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T};

opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fpd = {...                                        %when 2 outer rings are excluded!
    'FPD_Segmentation','OFF','FPD_Pixel',1,'FPD_MeanEff',0.95};%*0.8243243243};
%paper: 0.95 ?

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

opt_corr= {...
'RadiativeFlag','ON'};


% Tritium spectrum definition
myobj = TBD(...
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

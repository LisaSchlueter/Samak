function TritiumObject = ref_TBD_DScomparison(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration File for sensitivity studies 
% KATRIN nominal settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parser
p = inputParser;

% TBDDS 
p.addParameter('mnuSq_i',1,@(x)isfloat(x));
p.addParameter('Q_i',18575,@(x)isfloat(x)); % effective Endpoint

% WGTS
p.addParameter('WGTS_MolFrac_TT',0.95,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_DT',0.025,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_MolFrac_HT',0.025,@(x)isfloat(x)&& x>0);

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','WGTS100K','BlindingKNM1'}));
p.addParameter('HTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','SAENZ','WGTS100K','BlindingKNM1'}));
p.addParameter('TTFSD','BlindingKNM1',@(x)ismember(x,{'OFF','DOSS','SAENZ','WGTS100K','BlindingKNM1'}));

% KATRIN GENERAL SETTINGS
p.addParameter('TD','DScomparison',@(x) ischar(x));
p.addParameter('qU',linspace(18575-40,18575,10)',@(x)isfloat(x) && all(x>0));
p.addParameter('qUfrac',1/10.*ones(10,1)',@(x)isfloat(x));

% TD Reading Mode
p.addParameter('TDMode','DataTBD',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','FSD',@(x)ismember(x,{'OFF','numConv','matConv'}));

% Transmission Function
p.addParameter('KTFFlag','WGTSMACE',@(x)ismember(x,{'OFF','MACE','WGTSMACE'}));

p.parse(varargin{:});

% isotopologues concentrations
WGTS_MolFrac_HT          = p.Results.WGTS_MolFrac_HT;
WGTS_MolFrac_DT          = p.Results.WGTS_MolFrac_DT;
WGTS_MolFrac_TT          = p.Results.WGTS_MolFrac_TT;

% genereal
mnuSq_i           = p.Results.mnuSq_i;
Q_i               = p.Results.Q_i;
DopplerEffectFlag = p.Results.DopplerEffectFlag;
TD                = p.Results.TD;
TDMode            = p.Results.TDMode;
qU                       = p.Results.qU;
qUfrac                   = p.Results.qUfrac;

% final state distribution
DTFSD           = p.Results.DTFSD;
HTFSD           = p.Results.HTFSD;
TTFSD           = p.Results.TTFSD;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag; 

%%%%%%%%%%%%%%% END parser

opt_corr = {...
    'RadiativeFlag','ON',...
    'RecoilWmVmAFlag','OFF',...
    'FiniteExtChargeFlag','OFF',...
    'EEexchangeFlag','OFF',...    
    'ScreeningFlag','OFF',...
    'WintFiniteSizeFlag','OFF',...
    'RecoilCoulombFlag','OFF'};

opt_katrin = {...
    'TDMode',TDMode,...
    'TD',TD,...
    'mnuSq_i',mnuSq_i,...
    'Q_i',Q_i,...
    'qU',qU,...
    'qUfrac',qUfrac};

opt_wgts = {...
    'WGTS_MolFrac_TT',WGTS_MolFrac_TT,...
    'WGTS_MolFrac_DT',WGTS_MolFrac_DT,...
    'WGTS_MolFrac_HT',WGTS_MolFrac_HT};


opt_wgtsmace = {...
    'KTFFlag',KTFFlag};

opt_fsd= {...
    'TTFSD',TTFSD,...
    'DTFSD',DTFSD,...
    'HTFSD',HTFSD};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};


% Tritium spectrum definition
TritiumObject = TBD(...
    opt_corr{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_doppler{:},...
    opt_integration{:});
end

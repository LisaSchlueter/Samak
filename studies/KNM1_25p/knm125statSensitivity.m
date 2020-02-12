% KNM1 Neutrino Mass Sensitivity 
% 30 days of KATRIN
% 25% Column Density

% Data Taking Time
TimeDays = 54;
TimeSec = TimeDays*24*60*60;

% Pixel List for Un iform Fit
FPD_SegOffPixList = 1:119;
FPD_ROIlow        = 14;

% Column Density
WGTS_CD_MolPerCm2 = 1.1040e+17;
%WGTS_CD_MolPerCm2 = 5e+17;

% ELoss
ELossFlag         = 'KatrinD2';

% Magnetic Fields
MACE_Ba_T         = 6*1e-04;
WGTS_B_T          = 0.7*3.6;
MACE_Bmax_T       = 0.7*6;

% EndPoint
Q_i               = 18574;

% MTD
TD = 'KNM1_25p_40eV';

% FPD
FPD_ROIlow        = 14;
FPD_MeanEff       = 0.95;

% Analyis Range
range             = [30];
Scan              = 'OFF';
PlotFit           = 'ON';   %only used when Scan ON
SysBudget         = '06';    % 6=42 day, 07= 300 days, 09 = 900 days

% Background
BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',FPD_ROIlow,...
    'AnchorBkg6G',380e-3, ...
    'FPD_SegOffPixList',FPD_SegOffPixList);
                
%% Init Sensitivity Class
Arg = {...
    'TimeSec',TimeSec,...
    'FPD_SegOffPixList',FPD_SegOffPixList,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'Q_i',Q_i,...
    'range',range(1),...
    'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'Scan',Scan,...
    'PlotFit',PlotFit,...
    'SysBudget',SysBudget,...
    'TD',TD,...
    'mNuStop',1.5,...
    'FPD_MeanEff',FPD_MeanEff,...
    'ELossFlag',ELossFlag,...
    'FPD_ROIlow',FPD_ROIlow,...
    'BKG_RateSec',BKG_RateSec,...
    'chi2','chi2Stat'};

S = SensitivityStudy(Arg{:});

%%
S.NuMassScan

%% 
m90=908.3e-3;
sm2=m90^2/1.645;
m95=sqrt(2*sm2)*1000
m3sig=sqrt(3*sm2)*1000

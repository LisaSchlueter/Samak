function [CMObj, MultiCM, MultiCMFrac, MultiCMShape, MultiCMNorm]= ref_CovarianceMatrix_NominalKATRIN(varargin)
%This function does:
% Define systematic uncertainties
% Initialize Covariance Matrix Object
% Compute Covariance Matrix
%
% input: 
% neccessary:- TBD Object (ModelObj)
% optional:  - Sys uncertainties, sample size,...
% ouput:
% CovarianceMatrix Object
%
%            Lisa Schlüter  
%              TUM /MPP  
%           September 2018
%
% --------------------------- PARSER START ------------------------------------%
% Latest modifications:
%new 10.10.18 Lisa:
p = inputParser;
p.addParameter('SysBudget','01',@(x)ischar(x));
p.addParameter('ModelObj','',@(x)isa(x,'TBD'));
p.addParameter('PlotCM','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SysEffect','FSD'); %@(x)ischar(x) || @(x)isstruct(x)
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nTrials',1000,@(x)isfloat(x));

% default systematic uncertainties
p.addParameter('WGTS_TASR_RelErr',0,@(x)isfloat(x) && x>=0);
p.addParameter('FSDNorm_RelErr',0.01,@(x)isfloat(x) && x>=0);
p.addParameter('FSDShapeGS_RelErr',0.04,@(x)isfloat(x) && x>=0);
p.addParameter('FSDShapeES_RelErr',0.18,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Ba_T_RelErr',0.002,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Bmax_T_RelErr',0.002,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T_RelErr',0.002,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_CD_MolPerCm2_RelErr',0.001,@(x)isfloat(x) && x>=0);
p.addParameter('ISXsection_RelErr',0.001,@(x)isfloat(x) && x>=0);
p.addParameter('DTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('HTFlag','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

ModelObj                 = p.Results.ModelObj;
SysBudget                = p.Results.SysBudget;
PlotCM                   = p.Results.PlotCM; 
SysEffect                = p.Results.SysEffect;
RecomputeFlag            = p.Results.RecomputeFlag;
nTrials                  = p.Results.nTrials;
FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;
FSDShapeGS_RelErr        = p.Results.FSDShapeGS_RelErr;
FSDShapeES_RelErr        = p.Results.FSDShapeES_RelErr;
WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
MACE_Ba_T_RelErr         = p.Results.MACE_Ba_T_RelErr;
MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
ISXsection_RelErr        = p.Results.ISXsection_RelErr;
DTFlag                   = p.Results.DTFlag; 
TTFlag                   = p.Results.TTFlag; 
HTFlag                   = p.Results.HTFlag; 

if isempty(ModelObj)
    fprintf(2,'Input missing: ModelObj \n')
    return
end
% --------------------- END PARSER ---------------------------------%

switch SysBudget
    case '00'
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        MACE_Ba_T_RelErr = 0.002;
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.002;
        ISXsection_RelErr = 0.02;
        WGTS_TASR_RelErr = 0.001;
    case '01'
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        MACE_Ba_T_RelErr = 2*1e-06/ModelObj.MACE_Ba_T; %fixed to absolute value of 2*10^-6 T
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.001;
        ISXsection_RelErr = 0;
        WGTS_TASR_RelErr = 0;
    case '02'
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        MACE_Ba_T_RelErr = 2*1e-06/ModelObj.MACE_Ba_T; %fixed to absolute value of 2*10^-6 T
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.002;
        ISXsection_RelErr = 0.0046;
        WGTS_TASR_RelErr = 0;
   case '03' %for reproduction Marco Kleesieks values
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.001;
        ISXsection_RelErr = 0;
        WGTS_TASR_RelErr = 0;
    case '04'
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.002;
        ISXsection_RelErr = 0.0046;
        WGTS_TASR_RelErr = 0;
    case '05' % Closest to TDR
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        MACE_Ba_T_RelErr = 0.002;
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.0014;
        ISXsection_RelErr = 0.0014;
        WGTS_TASR_RelErr = 0;
    case '06' % Thierry - KATRIN 2019
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.01;
        WGTS_B_T_RelErr = 0.01;
        WGTS_CD_MolPerCm2_RelErr = 0.02;
        ISXsection_RelErr = 0.0046;
        WGTS_TASR_RelErr = 0.005;
    case '07' % Thierry - KATRIN 2021
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.005;
        WGTS_B_T_RelErr = 0.005;
        WGTS_CD_MolPerCm2_RelErr = 0.005;
        ISXsection_RelErr = 0.0046;
        WGTS_TASR_RelErr = 0.005;
    case '08' % Thierry - KATRIN 2023 - No E-loss
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.002;
        ISXsection_RelErr = 0.0046;
        WGTS_TASR_RelErr = 0.005;
    case '09' % KATRIN 2023 - No E-loss and smaller rhod sigma
        DTFlag='ON'; HTFlag = 'ON'; TTFlag = 'ON';
        FSDNorm_RelErr = 0.01;
        FSDShapeGS_RelErr = 0.04;
        FSDShapeES_RelErr = 0.18;
        %MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*3.9*1e-03+1.07*1e-06)/ModelObj.MACE_Ba_T;
        % NEW - WARNING - Thierry 23/2/2019
        MACE_Ba_T_RelErr = (ModelObj.MACE_Ba_T*1.87*1e-03 + 2.61*1e-06) /ModelObj.MACE_Ba_T;
        %
        MACE_Bmax_T_RelErr = 0.002;
        WGTS_B_T_RelErr = 0.002;
        WGTS_CD_MolPerCm2_RelErr = 0.001;
        ISXsection_RelErr = 0.0;
        WGTS_TASR_RelErr = 0;
end

if ischar(SysEffect)
    if strcmp(SysEffect,'RF') && ~strcmp(SysBudget,'08') && ~strcmp(SysBudget,'09')
        SysEffects = struct('RF_EL','ON','RF_BF','ON','RF_RX','ON');
    elseif  strcmp(SysEffect,'RF') % if SysBudget '08' or '09' NO ELOSS
         SysEffects = struct('RF_EL','OFF','RF_BF','ON','RF_RX','ON');
    elseif  strcmp(SysEffect,'RF_BFRX') %B-Fields + RhoD + inelSigma
        SysEffects = struct('RF_BF','ON','RF_RX','ON');
    elseif strcmp(SysEffect,'all') && ~strcmp(SysBudget,'08') && ~strcmp(SysBudget,'09')
        SysEffects = struct('RF_EL','ON','RF_BF','ON', 'RF_RX','ON', 'FSD','ON',...
            'TASR','OFF','TCoff_RAD','ON','TCoff_OTHER','ON');
    elseif strcmp(SysEffect,'all')
        SysEffects = struct('RF_EL','OFF','RF_BF','ON', 'RF_RX','ON', 'FSD','ON',...
            'TASR','OFF','TCoff_RAD','ON','TCoff_OTHER','ON');
    elseif strcmp(SysEffect,'TC')
        SysEffects = struct('TCoff_RAD','ON','TCoff_OTHER','ON');
    else
        SysEffects = struct(SysEffect, 'ON');
    end
else
    SysEffects = SysEffect;
end
% --------------------  Initialize Covariance Matrix----------------------------%
CMObj = CovarianceMatrix('StudyObject',ModelObj, 'nTrials',nTrials,...
    'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag,'SanityPlots','OFF',...
    'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
    'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,...
    'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
    'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
    'ISXsection_RelErr',ISXsection_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr,...
    'FSDNorm_RelErr',FSDNorm_RelErr,...
    'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr,...
    'DTFlag',DTFlag,'HTFlag',HTFlag,'TTFlag',TTFlag);
%Compute Covariance Matrix

CMObj.ComputeCM('PlotSaveCM',PlotCM,'SysBudget_Label',sprintf('SysBudget_%s',SysBudget));
if ischar(SysEffect) && strcmp(PlotCM,'ON')
    if contains(SysEffect,'RF')
        GetCMInfo(CMObj,SysBudget)
        CMObj.DisplayCMInfo('savename',sprintf('SysBudget_%s',SysBudget));
    end
end

%% 
% In case of Multipixel Fits
CMObj.StudyObject.ComputeTBDDS;
CMObj.StudyObject.ComputeTBDIS;
TBDIS_NoBKG = CMObj.StudyObject.TBDIS-CMObj.StudyObject.BKG_RateSec.*(CMObj.StudyObject.TimeSec.*CMObj.StudyObject.qUfrac); % column vector

CovMatFracCombiCell = repmat({CMObj.CovMatFrac},1,ModelObj.nPixels);
MultiCM = TBDIS_NoBKG.*blkdiag(CovMatFracCombiCell{:}).*TBDIS_NoBKG';

CovMatFracNormCell = repmat({CMObj.CovMatFracNorm},1,ModelObj.nPixels);
MultiCMNorm = TBDIS_NoBKG.*blkdiag(CovMatFracNormCell{:}).*TBDIS_NoBKG';
MultiCMFrac = blkdiag(CovMatFracCombiCell{:});

%CovMatFracShapeCell = repmat({CMObj.CovMatFracShape},1,ModelObj.nPixels);
%MultiCMShape = TBDIS_NoBKG.*blkdiag(CovMatFracShapeCell{:}).*TBDIS_NoBKG';

TBDIS_Sample = mvnrnd(TBDIS_NoBKG,MultiCM,5000);
            TBDIS_SumExpected = sum(TBDIS_NoBKG);
            TBDIS_SumSample = sum(TBDIS_Sample,2);
            TBDIS_SampleNorm = TBDIS_Sample.*(TBDIS_SumExpected./TBDIS_SumSample);
            MultiCMShape = cov(TBDIS_SampleNorm);
            
end
%% 
% p.addParameter('SysEffect','FSD',@(x)ischar(x));
% p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
% p.addParameter('nTrials',1000,@(x)isfloat(x));
% % default systematic uncertainties
% p.addParameter('WGTS_TASR_RelErr',0,@(x)isfloat(x) && x>=0);
% p.addParameter('FSDNorm_RelErr',0.01,@(x)isfloat(x) && x>=0);
% p.addParameter('FSDShapeGS_RelErr',0.04,@(x)isfloat(x) && x>=0);
% p.addParameter('FSDShapeES_RelErr',0.18,@(x)isfloat(x) && x>=0);
% p.addParameter('MACE_Bmax_T_RelErr',0.002,@(x)isfloat(x) && x>=0);
% p.addParameter('WGTS_B_T_RelErr',0.002,@(x)isfloat(x) && x>=0);
% p.addParameter('WGTS_CD_MolPerCm2_RelErr',0.001,@(x)isfloat(x) && x>=0);
% p.addParameter('ISXsection_RelErr',0,@(x)isfloat(x) && x>=0);
% p.addParameter('DTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
% p.addParameter('TTFlag','ON',@(x)ismember(x,{'ON','OFF'}));
% p.addParameter('HTFlag','ON',@(x)ismember(x,{'ON','OFF'}));

% ModelObj                 = p.Results.ModelObj;
% SysEffect                = p.Results.SysEffect;
% RecomputeFlag            = p.Results.RecomputeFlag;
% nTrials                  = p.Results.nTrials;
% MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
% WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
% WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
% ISXsection_RelErr        = p.Results.ISXsection_RelErr;
% WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
% FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;
% FSDShapeGS_RelErr        = p.Results.FSDShapeGS_RelErr;
% FSDShapeES_RelErr        = p.Results.FSDShapeES_RelErr;
% DTFlag                   = p.Results.DTFlag;
% TTFlag                   = p.Results.TTFlag;
% HTFlag                   = p.Results.HTFlag;
% PlotCM                   = p.Results.PlotCM; 


function CovarianceMatrixObject = ref_FT_CovarianceMatrix(obj,varargin)
% ------------------------------------------------------------------------%
% Get Covariance Object with default settings for First Tritium Campaign May/June 2018
% !! Please do not change value here !! Use RunAnalysis-ComputeCM(varargin)
%   
%  obj ----> RunAnalysis/MultiRunAnalysis Object
%
%                        L.Schlueter (06/2018)
% ------------------------------------------------------------------------%

defaultEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON',...  % Column Density, inel cross ection
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','ON',...
    'TCoff_OTHER','ON',...
    'DOPoff','OFF');

p = inputParser;
% settings
p.addParameter('SysEffects',defaultEffects,@(x)isstruct(x));
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nTrials',1000,@(x)isfloat(x));
% default systematic uncertainties
p.addParameter('WGTS_TASR_RelErr',0.023,@(x)isfloat(x) && x>=0);
p.addParameter('FSDNorm_RelErr',0.03,@(x)isfloat(x) && x>=0);
p.addParameter('FSDShape_RelErr',0.03,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Ba_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Bmax_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_CD_MolPerCm2_RelErr',0.1,@(x)isfloat(x) && x>=0);
p.addParameter('ISXsection_RelErr',0.02,@(x)isfloat(x) && x>=0);

p.parse(varargin{:});

SysEffects               = p.Results.SysEffects;
RecomputeFlag            = p.Results.RecomputeFlag;
nTrials                  = p.Results.nTrials;
MACE_Ba_T_RelErr         = p.Results.MACE_Ba_T_RelErr;
MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
ISXsection_RelErr        = p.Results.ISXsection_RelErr;
WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;
FSDShape_RelErr          = p.Results.FSDShape_RelErr;
%---------------------END CM parser----------------------------%

% Simplification for First Tritium (Campaign May/June 2018):
% Always use same (average) MTD for Covariance Matrix
% TD used in CovarianceMatrix.m only for Labeling!
% Runlist of all Tritium Runs
[ ~, ~ , FTRunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,char('40769.h5','40770.h5','40771.h5','40772.h5'));
if isempty(obj.RunNr)
    RunNr = 0;
elseif ~isempty(obj.RunNr)
    RunNr = obj.RunNr;
end

if ismember(RunNr,FTRunList) || all(ismember(obj.RunList,FTRunList))
    % set to average Stacked values (for CM only)
    obj.ModelObj.TD = 'RunStack_538_540_541_542_543_603_604_611_613_667_668_669_670_671_672_673_674_675_676_677_678_679_680_681_682_683_684_685_686_687_688_689_690_691_692_693ex2b';
    obj.ModelObj.WGTS_CD_MolPerCm2 = 4.453109415039991e+17;
end

% Create Covariance Matrix Object
CovarianceMatrixObject = CovarianceMatrix('StudyObject',obj.ModelObj, 'nTrials',nTrials,...
    'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag,'SanityPlots','OFF',...
    'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,...
    'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,...
    'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
    'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
    'ISXsection_RelErr',ISXsection_RelErr,...
    'WGTS_TASR_RelErr',WGTS_TASR_RelErr,...
    'FSDNorm_RelErr',FSDNorm_RelErr,'FSDShape_RelErr',FSDShape_RelErr);
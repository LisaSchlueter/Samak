function FitResult = MC_SensitivityStudy_FTKATRIN(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity Test for KATRIN First Tritium settings
% Compute simulated spectra (amount: nSamples)
% - Options: with Stat and/or Sys Fluctuations
% Fit these spectra
% - Options: -with or without systematics (Covariance Matrix, different Effects or all)
%            -some Parameters fixed (1=neutrino mass, 2=endpoint,3=Background, 4 = Normalization, 5&6= FSD Parameters) 
% Output = Fit Results, Simulated Spectra, Fit Spectra, ModelObject from Fit            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=inputParser;
%Simulation Model Parameter
p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SysFluct','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Q_i_Fit',18573.7,@(x)isfloat(x));
p.addParameter('fixPar','1 5 6',@(x)ischar(x));
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat','chi2CM'}));
p.addParameter('SysEffect','FSD',@(x)ischar(x)); %1 sys effect e.g. FSD
p.addParameter('nSamples',1000,@(x)isfloat(x));
p.addParameter('saveResults','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('plotFit','OFF',@(x)ismember(x,{'ON','OFF'}));
%Covariance Matrix settings
p.addParameter('DataDriven','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('WGTS_TASR_RelErr',0.005,@(x)isfloat(x) && x>=0);
p.addParameter('FSDNorm_RelErr',0.01,@(x)isfloat(x) && x>=0);
p.addParameter('FSDShapeGS_RelErr',0.04,@(x)isfloat(x) && x>=0);
p.addParameter('FSDShapeES_RelErr',0.18,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Ba_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('MACE_Bmax_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_CD_MolPerCm2_RelErr',0.05,@(x)isfloat(x) && x>=0);
p.addParameter('ISXsection_RelErr',0.02,@(x)isfloat(x) && x>=0);
p.parse(varargin{:});
DataDriven               = p.Results.DataDriven;
MACE_Ba_T_RelErr         = p.Results.MACE_Ba_T_RelErr;
MACE_Bmax_T_RelErr       = p.Results.MACE_Bmax_T_RelErr;
WGTS_B_T_RelErr          = p.Results.WGTS_B_T_RelErr;
WGTS_CD_MolPerCm2_RelErr = p.Results.WGTS_CD_MolPerCm2_RelErr;
ISXsection_RelErr        = p.Results.ISXsection_RelErr;
WGTS_TASR_RelErr         = p.Results.WGTS_TASR_RelErr;
FSDNorm_RelErr           = p.Results.FSDNorm_RelErr;
FSDShapeGS_RelErr        = p.Results.FSDShapeGS_RelErr;
FSDShapeES_RelErr        = p.Results.FSDShapeES_RelErr;
StatFluct                = p.Results.StatFluct;
SysFluct                 = p.Results.SysFluct;
Q_i_Fit                  = p.Results.Q_i_Fit;
fixPar                   = p.Results.fixPar;
exclDataStart            = p.Results.exclDataStart;
chi2                     = p.Results.chi2;
nSamples                 = p.Results.nSamples;
saveResults              = p.Results.saveResults;
plotFit                  = p.Results.plotFit;
SysEffect                = p.Results.SysEffect;
if strcmp(chi2,'chi2Stat')
    SysEffect = '';
end
%---------------------------------parser end---------------------------------------

%Init Model
MRA = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat',...%chi2 set later
                        'fixPar',fixPar,'exclDataStart',exclDataStart);
MRA.chi2 = chi2;
belowE0 = 18575-MRA.ModelObj.qU(exclDataStart);

%Covariance Matrix
if ~strcmp(chi2,'chi2Stat')
    if strcmp(SysEffect,'RF')
        SysEffects = struct('RF_EL','ON','RF_BF','ON','RF_RX','ON');
    elseif strcmp(SysEffect,'all')
       SysEffects = struct(...
                'RF_EL','ON',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON');
    elseif strcmp(SysEffect,'TC')
        SysEffects = struct('TCoff_RAD','ON','TCoff_OTHER','ON');
    else
        SysEffects = struct(SysEffect,'ON');
    end
    MRA.ComputeCM('SysEffects',SysEffects,'StackCM','OFF','DataDriven',DataDriven,...
        'InitNormFit','ON',...
        'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
        'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,'ISXsection_RelErr',ISXsection_RelErr,...
        'WGTS_TASR_RelErr',WGTS_TASR_RelErr,'FSDNorm_RelErr',FSDNorm_RelErr,...
        'FSDShapeGS_RelErr',FSDShapeGS_RelErr,'FSDShapeES_RelErr',FSDShapeES_RelErr);
elseif strcmp(chi2,'chi2Stat')
    MRA.InitModelObj_Norm_BKG;
    SysEffect = '';
end

%Compute MC Data set
switch StatFluct
    case 'ON'
        if strcmp(SysFluct,'OFF')
            TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.ModelObj.TBDIS',nSamples)'; % nSamples simulated integral spectra
        elseif strcmp(SysFluct,'ON')
            TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.FitCM+diag(MRA.ModelObj.TBDIS),nSamples)';
        end
    case 'OFF'
        if strcmp(SysFluct,'OFF')
            TBDIS_Sim = repmat(MRA.ModelObj.TBDIS,1,nSamples);
        elseif  strcmp(SysFluct,'ON')
            TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.FitCM,nSamples)';
        end
end
    
% Fit
MRA.ModelObj.Q_i = Q_i_Fit; %Set Endpoint value
FitResult = cell(nSamples,1);
TBDIS_Fit = zeros(MRA.ModelObj.nqU,nSamples);
for i=1:nSamples
 MRA.RunData.TBDIS = TBDIS_Sim(:,i);
 MRA.Fit;
 FitResult{i} = MRA.FitResult;
 TBDIS_Fit(:,i) = MRA.ModelObj.TBDIS;
 if strcmp(plotFit,'ON')
     MRA.PlotFit('LabelFlag','Simulation','Mode','Rate');
 end
end

% save Results to file
if strcmp(saveResults,'ON')
    if exist('./results/','dir')~=7 %in case results folder doesnt exist
        mkdir ./results/
    end
    save_name = sprintf('./results/MC_SensitivityStudy_FTKATRIN_%s_StatFluct%s_SysFluct%s_%s%s_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
        MRA.StackFileName,StatFluct,SysFluct,chi2,SysEffect,strrep(fixPar,' ',''),belowE0,Q_i_Fit*10,nSamples);
    save(save_name,'FitResult','TBDIS_Sim','TBDIS_Fit')
end
end
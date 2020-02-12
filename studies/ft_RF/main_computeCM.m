addpath(genpath('../../../Samak2.0'));
myEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON');  % RF Column Density, Cross Section
WGTS_CD_MolPerCm2_RelErr = 0.01*[1 2 5 10]'; % Fluctuation
WGTS_CD_MolPerCm2_local = 5e17.*[0.1 0.25 0.35 0.5 0.65 0.75 0.85 1]';

parfor i=1:numel(WGTS_CD_MolPerCm2_local)
 fprintf('Computing CM [CD var] %u out of %u', i,numel(WGTS_CD_MolPerCm2_local))
A = InitKatrin_ft_RF('WGTS_CD_MolPerCm2', WGTS_CD_MolPerCm2_local(i),'TD', 'FT-TL3');
C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,'SysEffect',myEffects,...
                      'WGTS_CD_MolPerCm2_RelErr',0.1,'RecomputeFlag','ON','SanityPlots','OFF');
C.ComputeCM_RF; %all ON
end

parfor j=1:numel(WGTS_CD_MolPerCm2_RelErr)
fprintf('Computing CM [CD var] %u out of %u', j,numel(WGTS_CD_MolPerCm2_RelErr))
 A = InitKatrin_ft_RF('WGTS_CD_MolPerCm2', 5e17 ,'TD', 'FT-TL4');
 C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,'SysEffect',myEffects,...
                      'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr(j),'RecomputeFlag','ON','SanityPlots','OFF');
 C.ComputeCM_RF; %all ON   
end
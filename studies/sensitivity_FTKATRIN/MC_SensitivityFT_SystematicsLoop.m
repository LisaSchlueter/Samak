% ------------------------------------------------------------------------------------------------
% Loop over MC_SensitivityFT with different Systematics
% Loop over qU-Scan ranges (200eV + 400eV)
% takes a while...to be run on server!
%----------------------------------------------------------     
% Goal:
% Obtain sensitivity on effective endpoint E0
% Method: Monte Carlo Simulation
% Simulate KATRIN nSamples times (asimov + add stat fluctuations)
% Options: 
% Fit with Statistic or use a covariance matrix , 3 Parameter free(E0,N,B)
% Plot E0 distribution
% Compute sensitivity to E0
% ------------------------------------------------------------------------------------------------
addpath(genpath('../../../Samak2.0'));
nSamples = 5000;
RunList = 'StackCD100all';
chi2 = 'chi2CM';
SysEffect  = {'TC','TASR','FSD','RF','all'};
fixPar = '1 5 6';
Q_i    = 18573.7;
exclDataStart_all = [7,9];
belowE0_all = [402,202];
FitResult = cell(nSamples,1);
RecomputeFlag = 'ON';

for j=1:2
    exclDataStart = exclDataStart_all(j);
    belowE0 = belowE0_all(j);
    % stat + sys
    parfor i=1:numel(SysEffect)
        save_name = sprintf('./results/MC_SensitivityStudy_FTKATRIN_%s_StatFluct%s_SysFluct%s_%s%s_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
            RunList,'ON','ON',chi2,SysEffect{i},strrep(fixPar,' ',''),belowE0,Q_i*10,nSamples);
        if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
            FitResult{i} = load(save_name)
        else
            FitResult{i} = MC_SensitivityStudy_FTKATRIN('saveResults','ON','StatFluct','ON','SysFluct','ON',...
                'chi2',chi2,'SysEffect',SysEffect{i},'nSamples',nSamples,'exclDataStart',exclDataStart);
        end
    end
    
    stat only
    save_name = sprintf('./results/MC_SensitivityStudy_FTKATRIN_%s_StatFluct%s_SysFluct%s_%s%s_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
        RunList,'ON','ON','chi2Stat','',strrep(fixPar,' ',''),belowE0,Q_i*10,nSamples);
    if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
        FitResultStat = load(save_name);
    else
        FitResultStat = MC_SensitivityStudy_FTKATRIN('saveResults','ON','StatFluct','ON','SysFluct','OFF',...
            'chi2','chi2Stat','nSamples',nSamples,'exclDataStart',exclDataStart);
    end
end
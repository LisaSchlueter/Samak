%-----------------------------------------------------------------------------------------
% Input: Choose 1 MTD (1 Ba, 1 range), Time intervals
% What this script does:
% 1) Calculate Nu-mass sensitivity using Scan Method, loop over systematic effects
% 2) Repeat 1) for different measurement times
% 3) Plot Result (work in progress...)
%-----------------------------------------------------------------------------------------
clear;
% Input:
SysBudget = '01'; %defines systematic budget in ref_CovarianceMatrix_...
MACE_Ba_T = 10*1e-04;
WGTS_B_T = 3.6*0.7;
range = 30;                                     % MTD range (eV below 18575)
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
TimeSec = (60*60*24*365).*[1/12, 2/12, 1/4, 1/2, 3/4, 1, 1.5, 2, 3, 3.5, 4, 5];
mnuSq_i_Fit = (0:20:200)*1e-03;                 % Neutrino mass scan range
SysEffect_all = {'TC','TASR','FSD','RF','all'};
RecomputeFlag = 'OFF';

% Init Gather results
chi2min = zeros(numel(SysEffect_all)+1,numel(mnuSq_i_Fit),numel(TimeSec));
mNu90   = zeros(numel(SysEffect_all)+1,numel(TimeSec)); %sensitivity of neutrino mass squared 90% C.L.
mNumin  = zeros(numel(SysEffect_all)+1,numel(TimeSec));
%% do stat only --> loop over time
SysEffect = '';
chi2 = 'chi2Stat';
save_file = sprintf('./results/%s_SensitivityStudy_NuMassScan_RunTime_%0feVrange_Ba%.0fG_%s%s.mat',SysBudget,range,MACE_Ba_T*1e4,chi2,SysEffect);
if exist(save_file,'file')==2 && strcmp(RecomputeFlag,'OFF')
    load(save_file);
else
   % progressbar('Computing...')
    parfor i=1:numel(TimeSec)
        %progressbar(i/numel(TimeSec));
        [~, chi2min(1,:,i), dof, mNu90(1,i),mNumin(1,i)] = SensitivityNominal_NuMassScan(...
            'chi2',chi2,'SysEffect',SysEffect,'TimeSec',TimeSec(i),'MACE_Ba_T',MACE_Ba_T,...
            'WGTS_B_T',WGTS_B_T,'TD',TD,'mnuSq_i_Fit',mnuSq_i_Fit,'plotFit','OFF','SysBuget',SysBudget,'CutScanRange','ON');
    end
    save(save_file,'TimeSec','mNu90','MACE_Ba_T','WGTS_B_T','mnuSq_i_Fit','TD','SysEffect','range','chi2min','chi2','dof');
end

%% do stat + systematics --> loop over sys + loop over time
chi2 = 'chi2CM';
for j=1:numel(SysEffects)
    SysEffect = SysEffect_all{i};
    save_file = sprintf('./results/SensitivityStudy_NuMassScan_RunTime_%0feVrange_Ba%.0fG_%s%s.mat',range,MACE_Ba_T*1e4,chi2,SysEffect);
    if exist(save_file,'file')==2 && strcmp(RecomputeFlag,'OFF')
        d = load(save_file);
    else
        sysIndex = j+1; %otherwise parfor doesnt work
        parfor ii=1:numel(TimeSec)
            [~, chi2min(sysIndex,:,ii), dof, mNu90(sysIndex,ii),mNumin(sysIndex,ii)] = SensitivityNominal_NuMassScan(...
                'chi2',chi2,'SysEffect',SysEffect,'TimeSec',TimeSec(ii),'MACE_Ba_T',MACE_Ba_T,...
                'WGTS_B_T',WGTS_B_T,'TD',TD,'mnuSq_i_Fit',mnuSq_i_Fit,'plotFit','OFF');
        end
        save(save_file,'TimeSec','mNu90','MACE_Ba_T','WGTS_B_T','mnuSq_i_Fit','TD','SysEffect','range','chi2min','dof');
    end
end
%% plot result
f70 = figure(70);
set(f70, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(TimeSec./(60*60*24*365),sqrt(mNu90)*1e3,'o--','LineWidth',3);
xlabel('Time (years)');
ylabel('\sigma_{m_{\nu}} 90% C.L. (meV)');
title(sprintf('Sensitivity Study (Scan Method) for KATRIN nominal settings\n stat + sys (%s) , %.0feV range ,  B_a = %.0fG , B_s = %.0fT',SysEffect,range,MACE_Ba_T*1e4,WGTS_B_T));
PrettyFigureFormat;
grid on;
set(gca,'FontSize',16);

save_name = sprintf('SensitivityNominal_NuMassScan_RunTime_%.0feVrange_Ba%.0fG',range,MACE_Ba_T*1e4);
export_fig(f70,['./plots/png/',save_name,'.png']);
savefig(f70,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f70,['./plots/pdf/',save_name,'.pdf']);
% Plot Sensitivity Curves for all systematic effects
% input: Scan range (below E0) and Ba
% output: plot
clear
range = 30;
MACE_Ba_T = 9*1e-04;
WGTS_B_T = 0.7*3.6;
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
SysBudget = '03';
save_name = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Systematics.mat',SysBudget,strrep(TD,'Sensitivity_',''),WGTS_B_T);
if exist(save_name,'file')==2
    load(save_name);
else
    fprintf('File doesnt exist! Run NuMassScan_SensitivityNominal_Systematics_Ba_Loop...\n')
    return
end
%% Plot of NuMassScanCurves
close
mNuSq_SysErr = sqrt(mNu90.^2-mNu90(1)^2);
stat_leg = ['Stat:                                                         \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(1)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(1))*1e3)];
TC_leg = ['Stat + Theoretical Corrections CM:         \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(2)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(2))*1e3)];  
TASR_leg = ['Stat + Tritium Activity Fluctuation CM:  \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(3)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(3))*1e3)];
FSD_leg = ['Stat + Final State Distribution CM:        \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(4)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(4))*1e3)];
RF_leg = ['Stat + Response Function CM:               \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(5)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(5))*1e3)];
All_leg = ['Stat + Combined CM:                             \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(6)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(6))*1e3)];

f33 = figure('Name','ScanSysOverview','Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
pScan = plot(mnuSq_i_Fit,chi2minScan,'LineWidth',3);
hold on;
ptitle = plot(mnuSq_i_Fit,NaN*zeros(numel(mnuSq_i_Fit),1),'Color',rgb('White'));
leg = legend([ptitle; pScan],'Sensitivity on neutrino mass (90% C.L.)',...
    stat_leg,TC_leg,TASR_leg,FSD_leg,RF_leg,All_leg,'Location','northwest');
legend boxoff
PrettyFigureFormat;
set(gca,'FontSize',18);
leg.FontSize = 14;
xlabel('m_{\nu^2} (eV^2)')
ylabel(sprintf('\\chi2 (%.0f dof)',dof(1)));
title(sprintf('Samak Fit to Simulation (asimov) \n - KATRIN %.0f years nominal settings - %.0feV range - Ba %.0fG, Bs %.2fT - BKG %.0f mcps',...
   TimeSec/(60*60*24*365),range,MACE_Ba_T*1e4,WGTS_B_T,GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T)*1e3));

save_name = sprintf('%sSysBudget_SensitivityNuMassScanOverview_%.0feVrange_Ba%0.fG_Bs%.2f',SysBudget,range,MACE_Ba_T*1e4,WGTS_B_T);
if ~exist('../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/','dir')
    mkdir ../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/
    mkdir ../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/
    mkdir ../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/
end
print(f33,['../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/',save_name,'.png'],'-dpng');
savefig(f33,['../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/',save_name,'.fig'],'compact');
publish_figurePDF(f33,['../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/',save_name,'.pdf']);






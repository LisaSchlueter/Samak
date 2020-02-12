clear;
MACE_Ba_T_all = [3:12]*1e-04;         % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
for b=1:numel(MACE_Ba_T_all)
MACE_Ba_T = MACE_Ba_T_all(b);
range = 60;                 % MTD energy range: 30,45,60eV below E0 (18575eV)
TimeSec = 3*365*24*60*60;
mnuSq_i_Fit = (0:20:200)*1e-03;
mySysEffects  = {'TC','TASR','FSD','RF','all'};
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
RecomputeFlag = 'OFF';

%Init: gather results
chi2min   = zeros(numel(mnuSq_i_Fit),numel(mySysEffects)+1); % chi2min distribution
mNu90     = zeros(numel(mySysEffects)+1,1);                  % sensitivity on mnu^2 (90% C.L.)     
mNumin    = zeros(numel(mySysEffects)+1,1);                  % mass with minimal chi2

save_name = sprintf('./results/SensitivityNominal_ResultsNuMassScan_MTD%s.mat',strrep(TD,'Sensitivity_',''));
if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
    load(save_name);
else    
[~, chi2min(:,1), dof, mNu90(1),mNumin(1)] = NuMassScan_SensitivityNominal(...
    'chi2','chi2Stat','TimeSec',TimeSec,'MACE_Ba_T',MACE_Ba_T,...
    'TD',TD,'mnuSq_i_Fit',mnuSq_i_Fit,'plotFit','OFF');

for i=1:numel(mySysEffects)
  [~, chi2min(:,i+1), ~, mNu90(i+1),mNumin(i+1)] = NuMassScan_SensitivityNominal(...
    'chi2','chi2CM','SysEffect',mySysEffects{i},...
    'TimeSec',TimeSec,'MACE_Ba_T',MACE_Ba_T,'TD',TD,...
    'mnuSq_i_Fit',mnuSq_i_Fit,'plotFit','OFF');  
end
save(save_name,'chi2min','mnuSq_i_Fit','mNu90', 'mNumin','TD','mySysEffects','dof');
end
% 
mNuSq_SysErr = sqrt(mNu90.^2-mNu90(1)^2);
%% Plot
close
stat_leg = ['Stat:                                                         \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(1)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(1))*1e3)];
TC_leg = ['Stat + Theoretical Corrections CM:         \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(2)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(2))*1e3)];  
TASR_leg = ['Stat + Tritium Activity Fluctuatoin CM:  \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(3)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(3))*1e3)];
FSD_leg = ['Stat + Final State Distribution CM:        \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(4)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(4))*1e3)];
RF_leg = ['Stat + Response Function CM:               \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(5)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(5))*1e3)];
All_leg = ['Stat + Combined CM:                             \sigma(m_{\nu^2}) = ',sprintf('%.3f eV² ,',mNu90(6)),' \sigma(m_{\nu}) = ', sprintf('%.0f meV',sqrt(mNu90(6))*1e3)];

f33 = figure(3);
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
pScan = plot(mnuSq_i_Fit,chi2min,'LineWidth',3);
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
title(sprintf('Samak Fit to Simulation (asimov) \n - KATRIN %.0f years nominal settings - %.0feV range - Ba %.0fG - BKG %.0f mcps',...
   TimeSec/(60*60*24*365),range,MACE_Ba_T*1e4,GetBackground(MACE_Ba_T)*1e3));

save_name = sprintf('SensitivityNuMassScanOverview_%.0feVrange_Ba%0.fG',range,MACE_Ba_T*1e4);
print(f33,['./plots/png/',save_name,'.png'],'-dpng');
savefig(f33,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f33,['./plots/pdf/',save_name,'.pdf']);
end

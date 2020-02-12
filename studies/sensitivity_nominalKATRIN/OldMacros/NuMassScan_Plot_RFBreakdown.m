% Output: numel(MACE_Ba_T) files 
clear;
SysBudget = '01'; 
WGTS_B_T = 3.6*0.7;
MACE_Ba_T = 7*1e-04;
range_all = [30,45,60];                        % MTD energy range: 30,45,60eV below E0 (18575eV)
d = cell(numel(range_all,1));
TimeSec = 3*(365*24*60*60);
Nranges = numel(range_all);
nSysRF = 5; 



%%
%Load File if possible
TD30 = sprintf('Sensitivity_%.0feV_Ba%.0fG',30,MACE_Ba_T*1e4);
TD60 = sprintf('Sensitivity_%.0feV_Ba%.0fG',60,MACE_Ba_T*1e4);
TD45 = sprintf('Sensitivity_%.0feV_Ba%.0fG',45,MACE_Ba_T*1e4);
save_name30 = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Systematics_RFcomponents.mat',SysBudget,strrep(TD30,'Sensitivity_',''),WGTS_B_T);
save_name60 = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Systematics_RFcomponents.mat',SysBudget,strrep(TD60,'Sensitivity_',''),WGTS_B_T);
save_name45 = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Systematics_RFcomponents.mat',SysBudget,strrep(TD45,'Sensitivity_',''),WGTS_B_T);
try
    d30 = load(save_name30);
    d45 = load(save_name45);
    d60 = load(save_name60);
catch
    fprintf('File doesnt exist! \n');
    return
end

mySysEffects = d30.mySysEffects;
s = GetSysBudet(SysBudget);
 %% plot
close all
f8 = figure(8);
set(f8, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
x = 1:1:numel(mySysEffects)+1;
mNu90Plot30 = sqrt(d30.mNu90)*1e3;
mNu90Plot30(2) = sqrt(d30.mNu90(4))*1e3; %swap energyloss and rhod
mNu90Plot30(4) = sqrt(d30.mNu90(2))*1e3;
mNu90Plot45 = sqrt(d45.mNu90)*1e3;
mNu90Plot45(2) = sqrt(d45.mNu90(4))*1e3; %swap energyloss and rhod
mNu90Plot45(4) = sqrt(d45.mNu90(2))*1e3;
mNu90Plot60 = sqrt(d60.mNu90)*1e3;
mNu90Plot60(2) = sqrt(d60.mNu90(4))*1e3; %swap energyloss and rhod
mNu90Plot60(4) = sqrt(d60.mNu90(2))*1e3;
plot(x,NaN*zeros(numel(x),1),'Color',rgb('White'));
hold on;
plot(x,mNu90Plot30,'--o','LineWidth',3,'Color',rgb('CornFlowerBlue'));
plot(x,mNu90Plot45,'--o','LineWidth',3,'Color',rgb('ForestGreen'));
plot(x,mNu90Plot60,'--o','LineWidth',3,'Color',rgb('IndianRed'));
xticks((1:1:numel(x)));
xticklabels({'Statistics only','Column Density + \sigma_{inel}', 'B-Fields','Energy Loss Function','All but Energy Loss','Response Function (all)'});
xtickangle(45);
ylim([200 350]);
yticks([200:10:350]);
ylabel('\sigma(m_\nu) 90% C.L. (meV)');
PrettyFigureFormat;
grid on;
legend('scan range','30 eV','45 eV','60 eV','Location','northwest');
title(sprintf('Response Function Uncertainty Breakdown - Nominal KATRIN %.0f years "optimized MTD" \n  B_a = %.0fG,  B_{T/max} = %.0f %%',...
    TimeSec/(365*24*60*60),MACE_Ba_T*1e4,WGTS_B_T/3.6*100));
legend boxoff
set(gca,'FontSize',14);
rf_leg = [sprintf('\\Delta \\rhod\\sigma = %.1f%%',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
     sprintf('\n\\DeltaB_a = %.0f \\muT, \\DeltaB_{T/max} = %.1f%%',1e6*s.MACE_Ba_T_Err,100*s.MACE_Bmax_T_RelErr),...
     sprintf('\nEnergy Loss Function (Uncertainties from Eur. Phys. J. D 10, 39–52)')];
 a=annotation('textbox',[0.47 0.27 1 0.1] , ...%upper corner: [0.155 0.755 1 0.1]
    'String', rf_leg, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left');
a.FontSize=11;a.FontWeight='bold';
save_name = sprintf('%s_SensitivityNominal_NuMassScan_RFBreakdown_Ba%.0fG_%.0fBFields',SysBudget,MACE_Ba_T*1e4,WGTS_B_T/3.6*100);
print(f8,['./plots/png/',save_name,'.png'],'-dpng');
savefig(f8,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f8,['./plots/pdf/',save_name,'.pdf']);

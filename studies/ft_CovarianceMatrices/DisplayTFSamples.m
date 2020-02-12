% Example of 1 Response Function Samples
RF_BF = 'ON'; RF_EL = 'ON'; RF_RX = 'ON';
nTrials = 1000;
NIS = 10;
TD = 'StackCD100allex2';
WGTS_CD_MolPerCm2 = 4.45982*1e17;
WGTS_CD_MolPerCm2_RelErr = 0.05;
ISXsection_RelErr = 0.005;
MACE_Ba_T = 6*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_B_T_RelErr = 0.02;
MACE_Ba_T_RelErr = 0.02;
if strcmp(RF_BF,'ON') && strcmp(RF_RX,'ON') && strcmp(RF_EL,'ON')
    rf_path = sprintf('../../inputs/CovMat/RF/LookupTables/RF/');
    rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr_Aseev.mat',...
        nTrials, NIS, TD, WGTS_CD_MolPerCm2,...
        WGTS_CD_MolPerCm2_RelErr,ISXsection_RelErr,...
        MACE_Bmax_T,WGTS_B_T,WGTS_B_T_RelErr,...
        floor(MACE_Ba_T*1e4),MACE_Ba_T_RelErr);
    rf_file = strcat(rf_path,rf_filename);
else
    rf_path = sprintf('../../inputs/CovMat/RF/LookupTables/PartVar/RFPartVar/');
    rf_filename = sprintf('RFInfo_%u-Trials_%u-NIS_BField-%s_RhoXsection-%s_Eloss-%s_%s_%g-molPercm2_%.2gerr_xsection%gerr_%.2fT-Bmax_%.2fT-Bs_BT%gerr_%.0fG-Ba_%.2gerr.mat',...
        nTrials, NIS,RF_BF,RF_RX,RF_EL,...
        TD,WGTS_CD_MolPerCm2,...
        WGTS_CD_MolPerCm2_RelErr,ISXsection_RelErr,...
        MACE_Bmax_T,WGTS_B_T,WGTS_B_T_RelErr,...
        floor(MACE_Ba_T*1e4),MACE_Ba_T_RelErr);
    rf_file = strcat(rf_path,rf_filename);
end
%rf_file = './DisplayInput/RFInfo_1000-Trials_11-NIS_StackCD100allex2_4.45982e+17-molPercm2_0.08err_xsection0.02err_4.20T-Bmax_2.52T-Bs_BT0.02err_6G-Ba_0.02err.mat';
load(rf_file); 

R = MultiRunAnalysis('RunList','StackCD100all');
%%
f53 = figure('Name','RF','Renderer','opengl');
set(f53, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
TF = zeros(R.ModelObj.nTe,1000);
for i=1:1000
TF(:,i) = R.ModelObj.ComputeMaceTF(R.ModelObj.Te,R.ModelObj.qU(15),'MACE_Ba_T',MACE_Ba_T_local(i),'WGTS_B_T',WGTS_B_T_local(i),...
    'MACE_Bmax_T',MACE_Bmax_T_local(i));
end

[l,a] = boundedline(R.ModelObj.Te-R.ModelObj.qU(15),mean(TF,2)*100,100*std(TF')');
l.LineWidth = 2; l.Color = rgb('CadetBlue');
a.FaceColor = rgb('CadetBlue'); a.FaceAlpha = 0.5;
xlim([-1 4]);
ylim([-5 105]);
PrettyFigureFormat;
set(gca,'FontSize',25);
xlabel('E_{kin} - retarding potential (V)');
ylabel('transmission probability (%)');
grid on;

save_name = sprintf('../ft_CovarianceMatrices/plots/ResponseFunction/TF1sigma');
print(f53,[save_name,'.png'],'-dpng');
publish_figurePDF(f53,[save_name,'.pdf']);

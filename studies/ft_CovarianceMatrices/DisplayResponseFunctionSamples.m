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

%% Reponse Function 1sigma error band
qUIndex = 5;
f5 = figure('Name','RF','Renderer','opengl');
set(f5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
[l, a] = boundedline(Te-qU(qUIndex),RFmean(qUIndex,:)*100,RFstd(qUIndex,:)*100);
l.LineWidth = 2; l.Color = rgb('CadetBlue');
a.FaceColor = rgb('CadetBlue'); a.FaceAlpha = 0.5;
PrettyFigureFormat;
xlim([-5 100]);
ylim([-3 100]);
legend(a,'1\sigma error band',...%sprintf('qU = %.0fV',qU(qUIndex)),
    'Location','northwest');
legend boxoff
xlabel('E_{kin} - retarding potential (V)');
ylabel('transmission probability (%)');
set(gca,'Fontsize',24);

save_name = sprintf('../ft_CovarianceMatrices/plots/ResponseFunction/RF_1sigma');
print(f5,[save_name,'.png'],'-dpng');
publish_figurePDF(f5,[save_name,'.pdf']);
%% varied parameter: example B source
WGTS_B_T_samples = WGTS_B_T.*(1+0.02.*randn(10000,1));
f52 = figure('Name','RF','Renderer','opengl');
set(f52, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6 ,0.6]);

h1 = histfit(WGTS_B_T_samples,20,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
PrettyFigureFormat;
xlabel('B_s (T)');
ylabel('samples');
yticks([0 500 1000 1500]);
legend(sprintf('<B_s> = %.2f T\n\\sigma(B_s) = %.1f %%',mean(WGTS_B_T_samples),std(WGTS_B_T_samples)/mean(WGTS_B_T_samples)*100));
legend boxoff;
set(gca,'Fontsize',25);

save_name = sprintf('../ft_CovarianceMatrices/plots/ResponseFunction/Bs_hist');
print(f52,[save_name,'.png'],'-dpng');
publish_figurePDF(f52,[save_name,'.pdf']);

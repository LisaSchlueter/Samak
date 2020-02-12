  % Example of 1 Response Function Samples
RF_BF = 'ON'; RF_EL = 'ON'; RF_RX = 'ON';
nTrials = 1000;
NIS = 11;
TD = 'StackCD100allex2';
WGTS_CD_MolPerCm2 = 4.45982*1e17;
WGTS_CD_MolPerCm2_RelErr = 0.08;
ISXsection_RelErr = 0.02;
MACE_Ba_T = 6*1e-04;
MACE_Bmax_T = 6*0.7;
WGTS_B_T    = 3.6*0.7;
WGTS_B_T_RelErr = 0.02;
MACE_Ba_T_RelErr = 0.02;

if strcmp(RF_BF,'ON') && strcmp(RF_RX,'ON')
    isp_file = sprintf('../../inputs/CovMat/RF/LookupTables/ISProb/ISProb-LookupTable_%u-Trials_%g-molPercm2-%.2gerr_xsection%gerr_%.0f-NIS_%.2f-Bmax_%.2f-Bs_BT%gerr.mat',...
        nTrials,WGTS_CD_MolPerCm2,WGTS_CD_MolPerCm2_RelErr,ISXsection_RelErr,NIS, MACE_Bmax_T,WGTS_B_T,WGTS_B_T_RelErr);
elseif strcmp(RF_BF,'ON') && strcmp(RF_RX,'OFF')
    isp_file = sprintf('../../inputs/CovMat/RF/LookupTables/PartVar/ISProbPartVar/ISProb-LookupTable_%u-Trials_%g-molPercm2_%.0f-NIS_%.2f-Bmax_%.2f-Bs_BT%gerr_BFFlag-%s.mat',...
        nTrials,WGTS_CD_MolPerCm2,NIS, MACE_Bmax_T, WGTS_B_T,WGTS_B_T_RelErr ,RF_BF);
elseif strcmp(RF_BF,'OFF') && strcmp(RF_RX,'ON')
    isp_file = sprintf('../../inputs/CovMat/RF/LookupTables/PartVar/ISProbPartVar/ISProb-LookupTable_%u-Trials_%g-molPercm2-%.2gerr_xsection%gerr_%.0f-NIS_%.2f-Bmax_%.2f-Bs_RXFlag-%s.mat',...
        nTrials,WGTS_CD_MolPerCm2,WGTS_CD_MolPerCm2_RelErr,ISXsection_RelErr,NIS, MACE_Bmax_T, WGTS_B_T, RF_RX);
elseif strcmp(RF_BF,'OFF') && strcmp(RF_RX,'OFF')
    isp_file = sprintf('../../inputs/WGTSMACE/WGTS_ISProb/IS_%g-molPercm2_%.0f-NIS_%.3f-Bmax_%.3f-Bs.mat',...
        WGTS_CD_MolPerCm2,NIS,MACE_Bmax_T, WGTS_B_T);
end
        
load(isp_file);
%% errorbars
close all
f122= figure('Renderer','opengl');
set(f122, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
b = bar(0:1:10,Pis_mean,'FaceColor',rgb('CadetBlue'),'LineStyle','none'); hold on;
e= errorbar(0:10,Pis_mean,std(Pis_m'),'.','Color',rgb('GoldenRod'),'LineWidth',2,'CapSize',15);
set(gca,'YScale','log');

ylabel('scattering probability (%)');
xlabel('number of scatterings');
PrettyFigureFormat;
set(gca,'FontSize',25);

save_name = sprintf('../ft_CovarianceMatrices/plots/ISProb/ISProb');
print(f122,[save_name,'.png'],'-dpng');
publish_figurePDF(f122,[save_name,'.pdf']);
%% correlation plot
f123= figure('Renderer','opengl');
set(f123, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.7]);
ISCovMat = cov(Pis_m');
corplot(ISCovMat)
PrettyFigureFormat;
yticklabels({'P_0','P_1','P_2','P_3','P_4','P_5','P_6','P_7','P_8','P_9','P_{10}'});
xticklabels({'P_0','P_1','P_2','P_3','P_4','P_5','P_6','P_7','P_8','P_9','P_{10}'});
c = colorbar; colormap;%(parula);
c.Label.String = 'correlation coefficient';
c.Label.FontSize = 25;
set(gca,'FontSize',25);
save_name = sprintf('../ft_CovarianceMatrices/plots/ISProb/ISProbCorPlot');
print(f123,[save_name,'.png'],'-dpng');
publish_figurePDF(f123,[save_name,'.pdf']);
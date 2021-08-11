

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savenameReg = sprintf('%sksn2_StatOverSyst_Twin.mat',savedir);
savenameFlat = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_LinFlatMTD.mat',savedir);
savenameIso = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_IsoStatMTD.mat',savedir);

dReg = importdata(savenameReg);
dFlat = importdata(savenameFlat);
dIso = importdata(savenameIso);

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.2,0.5]);
preg = plot(dReg.sin2t4_Stat.^2./dReg.sin2t4_Tot.^2,dReg.mNu4Sq);
hold on;
pflat = plot(dFlat.sin2T4_stat.^2./dFlat.sin2T4_cm.^2,dFlat.mNu4Sq);
piso = plot(dIso.sin2T4_stat.^2./dIso.sin2T4_cm.^2,dIso.mNu4Sq);
set(gca,'YScale','log');
PrettyFigureFormat;
leg = legend('ksn2','flat','iso stat','Location','southwest');
PrettyLegendFormat(leg);
xlim([0.5,1])
ylabel(sprintf('\\sigma_{stat.}^2 / \\sigma_{total}^2'))
xlabel(sprintf('{\\itm}_4^2 (eV^2)'));
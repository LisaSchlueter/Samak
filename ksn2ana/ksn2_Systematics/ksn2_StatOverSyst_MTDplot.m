

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savenameReg = sprintf('%sksn2_StatOverSyst_Twin.mat',savedir);
savenameRegTop3 = sprintf('%sksn2_StatOverSyst_Top3.mat',savedir);

savenameRegmBkgPT = sprintf('%sksn2_StatOverSyst_Sysm1_BkgPT.mat',savedir);

savenameFlat = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_LinFlatMTD.mat',savedir);
savenameIso = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_IsoStatMTD.mat',savedir);
% savenameReg2 = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_RegMTD.mat',savedir);
% savenameReg2_10mcps = sprintf('%sksn2_StatOverSyst_ref_KNM2_KATRIN_RegMTD_Bkg10mcps.mat',savedir);

dReg = importdata(savenameReg);
dRegTop3 = importdata(savenameRegTop3);

dRegmPT = importdata(savenameRegmBkgPT);
% dReg2 = importdata(savenameReg2);
% dReg2_10mcps = importdata(savenameReg2_10mcps);
dFlat = importdata(savenameFlat);
dIso = importdata(savenameIso);

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.2,0.5]);
pregall = plot(dReg.sin2t4_Stat.^2./dReg.sin2t4_Tot.^2,dReg.mNu4Sq,'LineWidth',2);
hold on;
pregTop3 = plot(dRegTop3.sin2T4_stat.^2./dRegTop3.sin2T4_cm.^2,dRegTop3.mNu4Sq,'LineWidth',2);
pregmPT = plot(dRegmPT.sin2T4_stat.^2./dRegmPT.sin2T4_cm.^2,dRegmPT.mNu4Sq,'LineWidth',2);
%preg = plot(dReg2.sin2T4_stat.^2./dReg2.sin2T4_cm.^2,dReg2.mNu4Sq,'LineWidth',2);
%preg10 = plot(dReg2_10mcps.sin2T4_stat.^2./dReg2_10mcps.sin2T4_cm.^2,dReg2_10mcps.mNu4Sq,'LineWidth',2);
pflat = plot(dFlat.sin2T4_stat.^2./dFlat.sin2T4_cm.^2,dFlat.mNu4Sq,'LineWidth',2);
piso = plot(dIso.sin2T4_stat.^2./dIso.sin2T4_cm.^2,dIso.mNu4Sq,'LineWidth',2);
set(gca,'YScale','log');
PrettyFigureFormat;
leg = legend([pregall,pregTop3,pregmPT,piso,pflat],...
    'ksn2 twin: all','ksn2 twin: top 3','ksn2 twin: all - PT',...
      ...  % 'ksn2',...%'ksn2 10 mcps',...
    'iso stat: top 3','flat: top 3','Location','southwest');
PrettyLegendFormat(leg);
xlim([0.5,1])
xlabel(sprintf('\\sigma_{stat.}^2 / \\sigma_{total}^2'))
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
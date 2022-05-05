% inconclusive, something wrong with isostat


savenameCommon = [getenv('SamakPath'),'ksn2ana/ksn2_MTD/results/Contour_E0NormBkg_KNM2_KATRIN_'];


dIso220 = importdata([savenameCommon,'IsoStatMTD.mat']);
dIso0 = importdata([savenameCommon,'IsoStatMTD_Bkg0mcps.mat']);
dReg = importdata([savenameCommon,'RegMTD_Backgrounds.mat']);


f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);

piso_reg = plot(dIso220.sin2T4_contour,dIso220.mNu4Sq_contour,'LineWidth',2,'Color',rgb('Red'),'LineStyle','-');
hold on;
piso_0 = plot(dIso0.sin2T4_contour,dIso0.mNu4Sq_contour,'LineWidth',2,'Color',rgb('Red'),'LineStyle','-.');
preg_reg = plot(dReg.sin2T4_contour_reg,dReg.mNu4Sq_contour_reg,'LineWidth',2,'Color',rgb('DodgerBlue'),'LineStyle','-');
preg_0 = plot(dReg.sin2T4_contour_0,dReg.mNu4Sq_contour_0,'LineWidth',2,'Color',rgb('DodgerBlue'),'LineStyle','-.');

set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);

leg = legend([preg_reg,piso_reg,preg_0,piso_0],...
    sprintf('KNM2 MTD, Bkg = 220 mcps'),...
     sprintf('IsoStat MTD, Bkg = 220 mcps'),...
      sprintf('KNM2 MTD, Bkg = 0 mcps'),...
       sprintf('IsoStat MTD, Bkg = 0 mcps'));
   PrettyLegendFormat(leg);
   leg.Position(2) = 0.6;
   
   xlim([6e-03, 0.5])
   ylim([1 1600]);
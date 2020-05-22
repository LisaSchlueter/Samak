% plot integral spectrum with and w/o sterile neutrino (exaggerated)
 RunAnaArg = {'RunList','KNM1',...
        'fixPar','E0 Norm Bkg',...
        'DataType','Twin',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'TwinBias_Q',18573.73,...
        'SysBudget',24,...
        'pullFlag',99,...
        'NonPoissonScaleFactor',1};
R = MultiRunAnalysis(RunAnaArg{:});
R.InitModelObj_Norm_BKG;
%% 

R.ModelObj.SetFitBiasSterile(0,0);
R.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
R.ModelObj.ComputeTBDDS; 
R.ModelObj.ComputeTBDIS; 
TimeFrac = R.ModelObj.qUfrac*R.ModelObj.TimeSec;
TBDIS_null = R.ModelObj.TBDIS./TimeFrac-R.ModelObj.BKG_RateSec;
TBDDS_null = R.ModelObj.TBDDS;

R.ModelObj.SetFitBiasSterile(50^2,0.3);
R.ModelObj.ComputeTBDDS; 
R.ModelObj.ComputeTBDIS; 
TBDIS_sterile = R.ModelObj.TBDIS./TimeFrac-R.ModelObj.BKG_RateSec;
TBDDS_sterile = R.ModelObj.TBDDS;

Te = R.ModelObj.Te; % energies electron
qU = R.ModelObj.qU; % retarding potentials
%% differential spectrum
NormFactor = TBDDS_null(end-7e2)./TBDDS_sterile(end-7e2);
GetFigure;
plot(Te-R.ModelObj.Q,TBDDS_null,'-','Color',rgb('DodgerBlue'),'LineWidth',2.5);
hold on;
plot(Te-R.ModelObj.Q,TBDDS_sterile.*NormFactor,'-.','Color',rgb('Orange'),'LineWidth',2.5);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('Electron energy - {\\itE}_0 (eV)'));
ylabel(sprintf('{\\itd}\\Gamma/{\\itdE} (cps \\cdot eV^{-1})'));
set(gca,'YScale','lin');
xlim([-90 0]);
ylim([0 25]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
legend('No sterile neutrino',sprintf('Sterile neutrino: {\\itm}_4 = 50 eV , {\\it sin}^2\\theta = 0.30'),...
    'EdgeColor',rgb('Silver'));
title('Differential spectrum','FontWeight','normal','FontSize',get(gca,'FontSize'));

plotdir = [getenv('SamakPath'),'ksn1ana/ksn1_systematics/plots/'];
MakeDir(plotdir);
plotname = sprintf('%sksn1_DiffSpecSterile.png',plotdir);
print(gcf,plotname,'-dpng','-r450');
%% integral spectrum
%NormFactor =TBDIS_null(end-7)./TBDIS_sterile(end-7);
GetFigure;
plot(qU-R.ModelObj.Q,TBDIS_null,'-','Color',rgb('DodgerBlue'),'LineWidth',2);
hold on;
plot(qU-R.ModelObj.Q,TBDIS_sterile.*NormFactor,'-.','Color',rgb('Orange'),'LineWidth',2);
hold on;
Diff = (TBDIS_sterile.*NormFactor-TBDIS_null);
Diff(Diff<0) = 0;
%plot(qU-R.ModelObj.Q,Diff,'-','Color',rgb('DimGray'),'LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('Retarding potential - {\\itE}_0 (eV)'))
ylabel(sprintf('Rate (cps)'));
set(gca,'YScale','lin');
xlim([-90 0]);
ylim([0 500])
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
legend('No sterile neutrino',sprintf('Sterile neutrino: {\\itm}_4 = 50 eV ,{\\it sin}^2\\theta = 0.30'),...
    'EdgeColor',rgb('Silver'));
title('Intgeral spectrum','FontWeight','normal','FontSize',get(gca,'FontSize'));
plotnameIS = sprintf('%sksn1_IntSpecSterile.png',plotdir);
print(gcf,plotnameIS,'-dpng','-r450');
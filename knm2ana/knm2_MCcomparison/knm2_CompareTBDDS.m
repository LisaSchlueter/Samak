% KNM2 Monte Carlo comparision
% Compare Differential spectrum
SanityPlot = 'ON';
ForeignModel = 'KaFit';
%% import Samak Diff Spec
savedir     = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
SamakStr    = [savedir,sprintf('knm2_ComputeTBDDS_NoFSD.mat')];
SamakStrFSD = [savedir,sprintf('knm2_ComputeTBDDS_FSD.mat')];

d  = importdata(SamakStr);
Te = d.Te;
TBDDS = d.TBDDS;
TBDDS = TBDDS./TBDDS(1);

dFSD  = importdata(SamakStrFSD);
TBDDS_FSD = dFSD.TBDDS;
TBDDS_FSD = TBDDS_FSD./TBDDS_FSD(1);
%% import Fitrium Diff Spec
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
FitriumStr = [savedir,'Fitrium_diffspec_NoFSD.dat'];
FitriumStrFSD = [savedir,'Fitrium_diffspec_FSD.dat'];

dF  = importdata(FitriumStr);
Te_F    = dF(:,1);
TBDDS_F = dF(:,2);
TBDDS_F = TBDDS_F./TBDDS_F(1);

dF_FSD  = importdata(FitriumStrFSD);
TBDDS_F_FSD = dF_FSD(:,2);
TBDDS_F_FSD = TBDDS_F_FSD./TBDDS_F_FSD(1);

%% import Kafit Diff Spec
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
KafitStr = [savedir,'KaFit_diffspec_NoFSD.dat'];
KafitStrFSD = [savedir,'KaFit_diffspec_FSD.dat'];

dK  = importdata(KafitStr);
Te_K    = dK(:,1);
TBDDS_K = dK(:,2);
TBDDS_K = TBDDS_K./TBDDS_K(1);

dK_FSD  = importdata(KafitStrFSD);
TBDDS_K_FSD = dK_FSD(:,2);
TBDDS_K_FSD = TBDDS_K_FSD./TBDDS_K_FSD(1);
%% %% plot difference
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 = plot(Te-18574,TBDDS-TBDDS_F,':','LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(Te-18574,TBDDS_FSD-TBDDS_F_FSD,'-','LineWidth',p1.LineWidth,'Color',rgb('Crimson'));
p3 = plot(Te-18574,TBDDS-TBDDS_K,'--','LineWidth',p1.LineWidth,'Color',rgb('Orange'));
p4 = plot(Te-18574,TBDDS_FSD-TBDDS_K_FSD,'-.','LineWidth',p1.LineWidth,'Color',rgb('Green'));
hold off;
leg = legend([p1,p3,p2,p4],'Samak - Fitrium  (No FSD)','Samak - KaFit     (No FSD)',...
    'Samak - Fitrium  (FSD)','Samak - KaFit     (FSD)');

leg.EdgeColor = rgb('Silver');
PrettyFigureFormat('FontSize',24)
xlabel('Energy - 18574 (eV)');
ylabel(sprintf('d\\Gamma/dE diff.'));
xlim([-91,10]);
plotdir = strrep(savedir,'results','plots');

savename = sprintf('%sTBDDS_Diff',plotdir);
%if strcmp(ForeignModel,'KaFit')
    ylim([0,5.2*1e-8]);
%end
export_fig(f2,[savename,'.pdf']);
print(f2,[savename,'.png'],'-dpng','-r300');

%% plot ratio
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plotY    = 1-TBDDS./y; plotY = plotY(plotY<0.2);
plotYFSD = 1-TBDDS_FSD./yFSD; plotYFSD = plotYFSD(plotYFSD<0.2);
plotX    = Te(plotY<0.2)-18574;
plotXFSD = Te(plotYFSD<0.2)-18574;
p3 = plot(Te-18574,1-TBDDS./TBDDS_K,'--','LineWidth',2.5,'Color',rgb('Orange'));
hold on;
p2 = plot(Te-18574,1-TBDDS_FSD./TBDDS_F_FSD,'-','LineWidth',p3.LineWidth,'Color',rgb('Crimson'));
p4 = plot(Te-18574,1-TBDDS_FSD./TBDDS_K_FSD,'-.','LineWidth',p3.LineWidth,'Color',rgb('Green'));
p1 = plot(Te-18574,1-TBDDS./TBDDS_F,':','LineWidth',p3.LineWidth,'Color',rgb('DodgerBlue'));
leg = legend([p1,p3,p2,p4],'Samak / Fitrium  (No FSD)','Samak / KaFit     (No FSD)',...
    'Samak / Fitrium  (FSD)','Samak / KaFit     (FSD)');
%leg = legend([p1,p2],['Samak/',ForeignModel,' (No FSD)'],['Samak/',ForeignModel,' (FSD)']);

leg.EdgeColor = rgb('Silver');
PrettyFigureFormat('FontSize',24)
xlabel('Energy - 18574 (eV)');
ylabel(sprintf('1 - rel. d\\Gamma/dE diff.'));
xlim([-91,10]);

leg.Location = 'southwest';
ylim([-4*1e-07,1*1e-7]);

savename = sprintf('%sTBDDS_Ratio',plotdir);
export_fig(f1,[savename,'.pdf']);
print(f1,[savename,'.png'],'-dpng','-r300');


% KNM2 response function comparison
SynchrotronFlag = 'ON';
AngTF = 'OFF';
%% Load RF Samak
RFConvBinStep = 0.01;
IntMode = 'Conv';
if strcmp(IntMode,'Conv')
    RFStepStr = sprintf('_RFbinStep%.2feV',RFConvBinStep);
elseif strcmp(IntMode,'Integral')
    RFStepStr = '';
end
savenameS = [getenv('SamakPath'),sprintf('knm2ana/knm2_MCcomparison/results/Response/Knm2_SamakRF%s_%s_Sync%s_ScatTF%s.mat',...
    RFStepStr,IntMode,SynchrotronFlag,AngTF)];

dSamak = importdata(savenameS);
qU = dSamak.qU;
TeS = dSamak.Te-qU; % energy
RFS = dSamak.RF; % transmission probability
fprintf('load file %s \n',savenameS)
%% Load KaFit
if strcmp(AngTF,'ON')
    savenameK = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Response/Knm2_KaFit_RF_detailed.dat'];
else
    savenameK = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Response/Knm2_KaFit_RF.dat'];
end
dKaFit = importdata(savenameK);
TeK = dKaFit(6:end-5,1); % energy
RFK = dKaFit(6:end-5,2); % transmission probability

%% Load Fitrium
savenameF = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Response/Knm2_Fitrium_RF.dat'];
dFitrium = importdata(savenameF);
TeF = dFitrium(:,1); % energy
RFF = dFitrium(:,2); % transmission probability

% %% plot ratio
% f2 = figure('Units','normalized','Position',[0,0.1,0.5,0.5]);
% l = plot(linspace(-5,90,100),ones(100,1),'-','LineWidth',2,'Color',rgb('Black'));
% hold on;
% pK = plot(TeS,RFS./RFK,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
% pF = plot(TeS,RFS./RFF,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
% hold on;
% pFconst = plot(TeS,RFS./RFFconst,':','LineWidth',2.5,'Color',rgb('IndianRed'));
% xlabel(sprintf('Energy - %.0f (eV)',qU));
% ylabel('Rel. probability diff.');
% PrettyFigureFormat('FontSize',22);
% leg = legend([pK,pF,pFconst],'Samak / KaFit','Samak / Fitrium',sprintf('Samak / Fitrium (const. \\sigma)'));
% leg.EdgeColor = rgb('Silver');
% xlim([min(TeS),max(TeS)]);
% 
plotdir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/plots/'];
% MakeDir(plotdir);
% savename = sprintf('%sRF_Ratio',plotdir);
% export_fig(f2,[savename,'.pdf']);
% print(f2,[savename,'.png'],'-dpng','-r300');

%% plot difference
f1 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.5]);
%hold on;
l = plot(linspace(-5,90,100),zeros(100,1),'-','LineWidth',2,'Color',rgb('Black'));
hold on;
pK = plot(TeS,RFS-RFK,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
pF = plot(TeS,RFS-RFF,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
pFK = plot(TeS,-RFF+RFK,':','LineWidth',2.5,'Color',rgb('IndianRed'));
%pFconst = plot(TeS,RFS-RFFconst,':','LineWidth',2.5,'Color',rgb('IndianRed'));
xlabel(sprintf('Energy - %.0f (eV)',qU));
ylabel('Probability diff.');
PrettyFigureFormat('FontSize',22);
leg = legend([pK,pF,pFK],'Samak - KaFit','Samak - Fitrium','KaFit - Fitrium');%,sprintf('Samak - Fitrium (const. \\sigma)'));
leg.EdgeColor = rgb('Silver');
leg.Location = 'northwest';
xlim([min(TeS),max(TeS)]);
%ylim([-3 5]*1e-04);
savename = sprintf('%sRF_Diff',plotdir);
%export_fig(f1,[savename,'.pdf']);
%print(f1,[savename,'.png'],'-dpng','-r300');


 %% ratio
 plotRatio = 'OFF';
 if strcmp(plotRatio,'ON') 
     f1 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.5]);
     l = plot(linspace(-5,90,100),zeros(100,1),'-','LineWidth',2,'Color',rgb('Black'));
     hold on;
     pK = plot(TeS,1-RFS./RFK,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
     pF = plot(TeS,1-RFS./RFF,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
     pFK = plot(TeS,1-RFK./RFF,':','LineWidth',2.5,'Color',rgb('IndianRed'));
     %pFconst = plot(TeS,1-RFS./RFFconst,':','LineWidth',2.5,'Color',rgb('IndianRed'));
     xlabel(sprintf('Energy - %.0f (eV)',qU));
     ylabel('1- Rel. probability diff.');
     PrettyFigureFormat('FontSize',22);
     leg = legend([pK,pF,pFK],'Samak / KaFit','Samak / Fitrium','KaFit / Fitrium');%,sprintf('Samak / Fitrium (const. \\sigma)'));
     leg.EdgeColor = rgb('Silver');
     leg.Location = 'northeast';
     xlim([min(TeS),max(TeS)]);
     ylim([-9 5]*1e-04);
     savename = sprintf('%sRF_Ratio',plotdir);
     export_fig(f1,[savename,'.pdf']);
     print(f1,[savename,'.png'],'-dpng','-r300');
 end
%% over lay
 
% f1 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.5]);
% l = plot(linspace(-5,90,100),zeros(100,1),'-','LineWidth',2,'Color',rgb('Black'));
% hold on;
% pS = plot(TeS,RFS,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));
% pK = plot(TeS,RFK,'-.','LineWidth',2.5,'Color',rgb('GoldenRod'));
% pF = plot(TeS,RFF,'--','LineWidth',2.5,'Color',rgb('ForestGreen'));
% pFconst = plot(TeS,RFFconst,':','LineWidth',2.5,'Color',rgb('IndianRed'));
% xlabel(sprintf('Energy - %.0f (eV)',qU));
% ylabel('Probability');
% PrettyFigureFormat('FontSize',22);
% leg = legend([pS,pK,pF,pFconst],'Samak','KaFit','Fitrium',sprintf('Fitrium (const. \\sigma)'));
% leg.EdgeColor = rgb('Silver');
% leg.Location = 'southeast';
% xlim([min(TeS),max(TeS)]);
% 
% savename = sprintf('%sRF_Overlay',plotdir);
% export_fig(f1,[savename,'.pdf']);
% print(f1,[savename,'.png'],'-dpng','-r300');

close all;
clear all;
savePlot = 'ON';

load('KNM1BestFit.mat');load('KNM2BestFit.mat');
LocalFontSize = 20;

fig1 = figure('Renderer','painters');
set(fig1, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

prlG = [81 126 102]/255;
prlB = [50 148 216]/255;
FitStyleArg1 = {'o','Color',rgb('Grey'),'LineWidth',1.0,'MarkerFaceColor',rgb('Grey'),'MarkerSize',4,'MarkerEdgeColor',rgb('Grey')};
FitStyleArg2 = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
err=sqrt(diag(B.FitCMShape));
err2=sqrt(diag(A.FitCMShape));

% Spectra  - First Subplot

s1=subplot(4,1,[1 2]);

% Plot

pfit = plot(A.RunData.qU(A.exclDataStart:end)-A.TwinBias_Q,...
    A.RunData.TBDIS(A.exclDataStart:end)./...
    (A.ModelObj.qUfrac(A.exclDataStart:end)*A.ModelObj.TimeSec),...
    'DisplayName','KNM1','color',prlG,'LineWidth',3,'LineStyle','-');
hold on

pdata = errorbar(A.RunData.qU(A.exclDataStart:end)-A.TwinBias_Q,...
    A.RunData.TBDIS(A.exclDataStart:end)./...
    (A.ModelObj.qUfrac(A.exclDataStart:end)*A.ModelObj.TimeSec),...
    (A.RunData.TBDISE(A.exclDataStart:end))./(A.ModelObj.qUfrac(A.exclDataStart:end)*A.ModelObj.TimeSec)*50,FitStyleArg1{:},'CapSize',0);

pfit2= plot(B.RunData.qU(B.exclDataStart:end)-B.TwinBias_Q,...
    B.RunData.TBDIS(B.exclDataStart:end)./...
    (B.ModelObj.qUfrac(B.exclDataStart:end)*B.ModelObj.TimeSec),...
    'DisplayName','KNM2','color',prlB,'LineWidth',3,'LineStyle','-');

pdata2= errorbar(B.RunData.qU(B.exclDataStart:end)-B.TwinBias_Q,...
    B.RunData.TBDIS(B.exclDataStart:end)./...
    (B.ModelObj.qUfrac(B.exclDataStart:end)*B.ModelObj.TimeSec),...
    (B.RunData.TBDISE(B.exclDataStart:end))./(B.ModelObj.qUfrac(B.exclDataStart:end)*B.ModelObj.TimeSec)*50,FitStyleArg2{:},'CapSize',0);

yl1 = ylabel('Count rate (cps)');
%annotation('rectangle',[0.18 0.888 0.05 0.003],'FaceColor',rgb('Grey'));
legend([pfit,pfit2],{'Spectrum KRN1\newline with 1\sigma error bars \times 50','Spectrum KRN2\newline with 1\sigma error bars \times 50'},'Location','northeast','box','off');
lgd=legend;
lgA.FontSize = LocalFontSize-2;

PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
set(gca, 'YScale', 'log');

s2=subplot(4,1,3);
boundedline(B.RunData.qU(B.exclDataStart:end)-B.TwinBias_Q,zeros(1,numel(B.RunData.qU(B.exclDataStart:end))),ones(1,numel(B.RunData.qU(B.exclDataStart:end))));
hold on;
res1=errorbar(B.RunData.qU(B.exclDataStart:end)-B.TwinBias_Q,...
    (B.RunData.TBDIS(B.exclDataStart:end)-B.ModelObj.TBDIS(B.exclDataStart:end))./err(B.exclDataStart:end),...
    zeros(1,numel(B.RunData.qU(B.exclDataStart:end))),FitStyleArg2{:});
res2=errorbar(A.RunData.qU(A.exclDataStart:end)-A.TwinBias_Q,...
    (A.RunData.TBDIS(A.exclDataStart:end)-A.ModelObj.TBDIS(A.exclDataStart:end))./err2(A.exclDataStart:end),...
    zeros(1,numel(A.RunData.qU(A.exclDataStart:end))),FitStyleArg1{:});
yl2=ylabel('Residuals (\sigma)');
legend([res1 res2],{'KRN1','KRN2'},'Location','southeast');
ylim([-2.5 2.5]);
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
hl.Position(2) = 0.333;
ax2 = gca;

% MTD      -  Third Subplot
s3=subplot(4,1,4);
times1 = A.ModelObj.qUfrac*A.ModelObj.TimeSec;
times2 = B.ModelObj.qUfrac*B.ModelObj.TimeSec;
bar(A.RunData.qU-A.TwinBias_Q,times1./(60*60*24),0.5,...
'FaceColor',prlG,'EdgeColor',prlG);
hold on;
bar(B.RunData.qU-B.TwinBias_Q,times2./(60*60*24),0.5,...
'FaceColor',prlB,'EdgeColor',prlB,'FaceAlpha',0.5,'EdgeAlpha',0.5);
legend({'KRN1','KRN2'});
xlabel('Retarding energy - 18575 (eV)');
ylh = ylabel('Time (d)');
ylh.Position(1) = ax2.YLabel.Position(1);%-range-6.8;
yl1.Position(1) = ax2.YLabel.Position(1);
yl2.Position(1) = ax2.YLabel.Position(1);
yl1.Position(2) = ax2.YLabel.Position(2)+4;
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);     

linkaxes([s1,s2,s3],'x');
xlim([-41 140]);
hold off;

fig2=figure('Renderer','painters');
set(fig2, 'Units', 'normalized', 'Position', [0.001, 0.001, 0.45, 0.37]);
Y=[2e14 1e13 7e11 1.12e11 1.29e11];
Yerr=[0 0 1.4e11 3.7e10 3.3e10];
FitStyleArg = {'o','Color',rgb('Red'),'LineWidth',1.0,'MarkerFaceColor',rgb('Red'),'MarkerSize',8,'MarkerEdgeColor',rgb('Red'),'CapSize',0};
errorbar([1 2 3 4 5],Y,Yerr,FitStyleArg{:});
xticks([1 2 3 4 5]);
xticklabels({'Los Alamos','Troitsk','KRN1','KRN2','Combined'});
ylabel('\eta','FontSize',12);
ax=gca;
ax.YScale='log';
xlim([0.5 5.5]);
ylim([1e5 1e15]);
grid on;
b=annotation('rectangle',[0.13 0.106 0.775 0.085],'FaceColor','blue','FaceAlpha',.2);
text(2,3e6,{'Limit from Pauli blocking'},'FontSize',20,'Color','black');
c=annotation('line',[0.125 0.9],[0.54 0.54]);
c.LineWidth=2;
c.LineStyle='--';
text(2.25,4e9,{'KATRIN sensitivity'},'FontSize',20,'Color','black');
PrettyFigureFormat;

fig3=figure('Renderer','Painters');
set(fig3, 'Units', 'normalized', 'Position', [0.001, 0.001, 0.45, 0.4]);
load('EtaFitResult_KNM1_AllParams_mnuSq0_Nfit1');D1=D;M1=M;load('EtaFitResult_KNM2_Prompt_AllParams_mnuSq0_Nfit1');
times = M1.ModelObj.qUfrac*M1.ModelObj.TimeSec;
qU    = M1.ModelObj.qU; qU    = qU-M1.TwinBias_Q; % Energy axis
IS  = D1.ModelObj.TBDIS; 
YIs = IS./times;
YI = M1.ModelObj.TBDIS; 
YI = YI./times;
DIS = D1.RunData.TBDIS;
DIS = DIS./times;
err  = (diag(sqrt(D1.FitCMShape)));
err  = err./times;
err  = err./YI;
s1=subplot(2,1,1);

% Ratio
RSP  = (YI./YIs);
RSPd = RSP;

% Plot
hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
hold on;
hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
hr3 = errorbar(qU,(DIS./YIs),err,FitStyleArg2{:},'CapSize',0);
yl2=ylabel(sprintf('Ratio'));
ax=gca;
katrinsim   = sprintf('\\eta=0');
sterilemod  = sprintf('Best fit: \\eta=%.2g',M1.ModelObj.eta);
hl1=legend([hr2 hr1],{katrinsim,sterilemod},'Location','southeast','box','off');
hl1.NumColumns=2;
hl1.FontSize = LocalFontSize-2;

xlim([-45 40]);
%ylim([min((DIS./YI-1)./err) max((DIS./YI-1)./err)]);
title('KRN1','FontWeight','normal');
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
hl.Position(2) = 0.333;
ax2 = gca;
xticks([]);
times = M.ModelObj.qUfrac*M.ModelObj.TimeSec;
qU    = M.ModelObj.qU; qU    = qU-M.TwinBias_Q; % Energy axis
qU    = M.ModelObj.qU; qU    = qU-M.TwinBias_Q; % Energy axis
IS  = D.ModelObj.TBDIS; 
YIs = IS./times;
YI = M.ModelObj.TBDIS; 
YI = YI./times;
DIS = D.RunData.TBDIS;
DIS = DIS./times;
err  = (diag(sqrt(D.FitCMShape)));
err  = err./times;
err  = err./YI;
s2=subplot(2,1,2);

% Ratio
RSP  = (YI./YIs);
RSPd = RSP;

% Plot
hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
hold on;
hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
hr3 = errorbar(qU,(DIS./YIs),err,FitStyleArg2{:},'CapSize',0);
%yl=ylabel(sprintf('2^{nd}'));
katrinsim   = sprintf('\\eta=0');
sterilemod  = sprintf('Best fit: \\eta=%.2g',M.ModelObj.eta);
hl=legend([hr2 hr1],{katrinsim,sterilemod},'Location','best','box','off');
hl.NumColumns=2;
hl.FontSize = LocalFontSize-2;
title('KRN2','FontWeight','normal');
xlim([-45 40]);
ylim([0.97 inf]);
xlabel('Retarding energy - 18575 (eV)');
PRLFormat;
set(gca,'FontSize',LocalFontSize);
set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
ax2 = gca;
yl2.Position(2) = ax.YLabel.Position(2)-0.04;

fig4=figure('Renderer','painters');
set(fig4, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.6]);
hold on;
load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq0_SystON_range40_KNM1_2D_Real.mat');
[~,u]=contour(mnuScanPoints(5:end),etaScanPoints,(Chi2(5:end,:)-GlobalChi2Min).',[0 9.21],'LineStyle',':');
u.LineWidth=4;
u.Color=prlG;
load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq0_SystON_range40_KNM2_Prompt_2D_Real.mat');
[~,v]=contour(mnuScanPoints(5:end),etaScanPoints,(Chi2(5:end,:)-GlobalChi2Min).',[0 9.21],'LineStyle','--');
v.LineWidth=4;
v.Color=prlB;
load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq0_SystON_range40_KNM2_Prompt_2D_Real_CombiFit.mat');
Chi2_2 = Chi2;
load('./RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq0_SystON_range40_KNM1_2D_Real_CombiFit.mat');
Chi2Combi = Chi2 + Chi2_2;
ZeroPoint = Chi2Combi(mnuScanPoints==min(abs(mnuScanPoints)),etaScanPoints==min(abs(etaScanPoints)));
[~,w]=contour(mnuScanPoints(1:end),etaScanPoints,(Chi2Combi(1:end,:)-ZeroPoint).',[0 9.21],'LineStyle','-');
w.LineWidth=4;
w.Color='red';
o=plot([0 1],[1.8e10 1.8e10],'LineWidth',4,'LineStyle','-.','Color','black');
xlim([0 1]);
ylim([1e10 inf]);
ax=gca;ax.YScale = 'log';
xlabel('m_\nu^2 (eV^2)','FontSize',12);
ylabel('\eta','FontSize',12);
l=plot([0 0],'Color',prlG,'LineWidth',4,'LineStyle',':');m=plot([0 0],'Color',prlB,'LineWidth',4,'LineStyle','--');n=plot([0 0],'Color','red','LineWidth',4,'LineStyle','-');
leg=legend([l m n o],{'KRN1','KRN2','Combination','KATRIN sensitivity limit'},'Location','none','box','off');
leg.Position=[0.47 0.3 0.5 0.2];
grid on;
PRLFormat;
hold off;

fig5=figure('Renderer','painters');
set(fig5, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.5]);
A=ref_RelicNuBkg_KNM1('TTFSD','OFF','HTFSD','OFF','DTFSD','OFF','mnuSq_i',0.5);
B=ref_RelicNuBkg_KNM1('mnuSq_i',0.5);
A.ComputeTBDDS;B.ComputeTBDDS;
a=semilogy(A.Te-A.Q+1.7,A.TBDDS+1e-50,'LineWidth',2);
hold on;
d=semilogy(B.Te-B.Q+1.7,B.TBDDS_R,'LineWidth',2);
c=semilogy(B.Te-B.Q+1.7,B.TBDDS-B.TBDDS_R+1e-50,'LineWidth',2);
b=semilogy(A.Te-A.Q+1.7,A.TBDDS_R,'LineWidth',2);
ylabel('Counts per 0.1 eV');
xlabel('E-E_0 (eV)');
legend([a b c d],'\beta decay, \zeta_{j}(V_{j})=\delta(0)','C\nuB, \zeta_{j}(V_{j})=\delta(0)','\beta decay, standard FSD','C\nuB, standard FSD');
legend boxoff;
xlim([-6+1.7 4+1.7]);
ylim([1e-24 1]);
x1=[0.651 0.651];
y1=[0.6 0.15];
a=annotation('line',x1,y1);
a.LineStyle='--';
a.LineWidth=2;
a.Color=[0.5 0.5 0.5];
x2=[0.517 0.517];
y2=[0.6 0.15];
b=annotation('line',x2,y2);
b.LineStyle='--';
b.LineWidth=2;
b.Color=[0.5 0.5 0.5];
x3=[0.613 0.651];
y3=[0.38 0.38];
c=annotation('textarrow',x3,y3,'String','2m_\nu');
c.FontSize=16;
x4=[0.56 0.538];
y4=[0.38 0.38];
d=annotation('arrow',x4,y4);
x5=[0.62 0.651];
y5=[0.55 0.55];
e=annotation('textarrow',x5,y5,'String','\langle E_{GS}\rangle');
e.FontSize=16;
x6=[0.547 0.517];
y6=[0.55 0.55];
f=annotation('arrow',x6,y6);
PrettyFigureFormat;
hold off;

A=ref_RelicNuBkg_DesignReport('ToggleES','ON','eta_i',1,'mnuSq_i',0.5);
B=ref_RelicNuBkg_DesignReport('ToggleES','OFF','eta_i',1,'mnuSq_i',0.5);
A.ComputeTBDDS;B.ComputeTBDDS;A.ComputeTBDIS;B.ComputeTBDIS;
fig6=figure('Renderer','painters');
set(fig6, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.45]);
semilogy(A.Te-A.Q+1.7,A.TBDDS_R,'LineWidth',2);
hold on;
plot(B.Te-B.Q+1.7,B.TBDDS_R,'LineWidth',2);
xlabel('E-E_0 (eV)');
ylabel('Counts per 0.1 eV');
xlim([-30-1.7 5]);
ylim([5e-16 8e-14]);
legend('Neutrino capture full differential spectrum','Gaussian fit of ground state (used in the analysis)','box','off','location','northwest');
PrettyFigureFormat;
hold off;

if strcmp(savePlot,'ON')
    SaveDir = [getenv('SamakPath'),sprintf('../../home/fkellere/Dokumente/KRN1 paper/Figures/')];
    MakeDir(SaveDir);
    SaveName1  = 'DataSets.pdf';
    SaveName2  = 'Overview.pdf';
    SaveName3  = 'BestFits.pdf';
    SaveName4  = 'Contours.png';
    SaveName5  = 'FSDeffects.pdf';
    SaveName6  = 'CnuB_DiffSpec.pdf';
    export_fig(fig1, [SaveDir,SaveName1]);
    export_fig(fig2, [SaveDir,SaveName2]);
    export_fig(fig3, [SaveDir,SaveName3]);
    export_fig(fig4, [SaveDir,SaveName4],'-m10');
    export_fig(fig5, [SaveDir,SaveName5]);
    export_fig(fig6, [SaveDir,SaveName6]);
    close all;
end
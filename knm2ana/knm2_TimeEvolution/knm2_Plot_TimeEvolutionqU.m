% plot mean retarding potential for knm1
% for PhD thesis
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedir,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
d = importdata(savename);

FontSize = 22;
close all

GetFigure;
s1 = subplot(4,1,[1:3]);
plot(linspace(18400,18640,10),zeros(10,1),'-','LineWidth',2,'Color',rgb('Silver'));
hold on;
for i=1:d.A.ModelObj.nqU
dscatter(1e-10.*rand(d.A.nRuns,1)+repmat(d.A.RunData.qU(i),1,d.A.nRuns)',...
    1e3*(d.A.SingleRunData.qU(i,:)-mean(d.A.SingleRunData.qU(i,:)))',...
    'Msize',10,'Marker','o');
end

ax = gca;
xticklabels(''); 
ylabel('Scan-wise variation (mV)');%sprintf('U - \\langle U\\rangle (mV)'));
yticks(1e3*[-0.1:0.1:0.1]);
PrettyFigureFormat('FontSize',22);
c = colorbar;
c.Label.String = 'Rel. scan density'; c.Label.FontSize = ax.XLabel.FontSize;
c.Position(3) = 0.01;c.Position(1) = 0.91;
colormap(flipud(colormap('winter')));


qU = d.A.SingleRunData.qU;

s2 = subplot(4,1,4);
bar(d.A.RunData.qU,1e3*std(qU-mean(qU,2),0,2),'FaceColor',rgb('Orange'),'EdgeColor',rgb('Orange'));
hold on;
pm = plot(linspace(18400,18640,10),1e3*mean(std(qU-mean(qU,2),0,2)).*ones(10,1),':','LineWidth',2,'Color',rgb('FireBrick'));
PrettyFigureFormat('FontSize',22);


ax2 = gca;
ax2.XAxis.Exponent = 0;
xlabel('Mean retarding energy (eV)');
ylabel(sprintf('\\sigma (mV)'));

linkaxes([s1,s2],'x');
xlim([18477 18625]);
ylim([0 60])
leg = legend(pm,sprintf('\\langle\\sigma\\rangle = %.0f mV',1e3*mean(std(qU-mean(qU,2),0,2))));
legend box off
ax2.Position(3) = ax.Position(3);
ax2.Position(4) = 0.19;
ax2.YLabel.Position(1) = ax.YLabel.Position(1);

pltdir = [getenv('SamakPath'),'knm2ana/knm2_TimeEvolution/plots/'];
MakeDir(pltdir);
pltfile = sprintf('%sknm2_qU_TimeEvolution.pdf',pltdir);
export_fig(pltfile);

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
plot(d.A.SingleRunData.StartTimeStamp,zeros(361,1),'-','LineWidth',2,'Color',rgb('Silver'));
hold on;
plot(d.A.SingleRunData.StartTimeStamp,-100.*ones(361,1),':','LineWidth',2,'Color',rgb('Silver'));
plot(d.A.SingleRunData.StartTimeStamp,100.*ones(361,1),':','LineWidth',2,'Color',rgb('Silver'));

p1 = plot(d.A.SingleRunData.StartTimeStamp,1e3.*(qU(32,:)-d.A.RunData.qU(32)),'.','MarkerSize',10,'Color',rgb('Crimson'));
%p2 = plot(d.A.SingleRunData.StartTimeStamp,1e3.*(qU(8,:)-d.A.RunData.qU(8)),'.','MarkerSize',10,'Color',rgb('DodgerBlue'));

MaxD = zeros(38,1);
 for i=1:38
     MaxD(i) = max(abs(1e3.*(qU(i,:)-d.A.RunData.qU(i))));
%     p1 = plot(d.A.SingleRunData.StartTimeStamp,1e3.*(qU(i,:)-d.A.RunData.qU(i)),'.','MarkerSize',10);
 end
ax = gca;
ax.YAxis.Exponent = 0;
ylim([-135 260])
ylabel('Scan-wise variation (mV)');
xlabel(sprintf('Date in %s',datestr(d.A.SingleRunData.StartTimeStamp(1),'yyyy')));
a = xticks;
xticklabels(datestr(a,'mmm, dd'))
leg = legend([p1],sprintf('\\langleqU\\rangle = %.1f eV',d.A.RunData.qU(32)));%,...
   % sprintf('\\langleqU\\rangle = %.1f eV',d.A.RunData.qU(8)));
PrettyLegendFormat(leg);
PrettyFigureFormat('FontSize',22);
leg.FontSize = get(gca,'FontSize')+2;

pltdir = [getenv('SamakPath'),'knm2ana/knm2_TimeEvolution/plots/'];
pltfile = sprintf('%sknm2_qU_TimeEvolution_Outlier.pdf',pltdir);
export_fig(pltfile);

%% outlier
fprintf('Outlier 1: <qU> = %.1f eV,  all std = %.0f meV , after %s std = %.0f meV \n',...
    d.A.RunData.qU(32),1e3.*std(d.A.SingleRunData.qU(32,:)),...
    d.A.SingleRunData.StartTimeStamp(46),1e3.*std(d.A.SingleRunData.qU(32,46:end)));

% fprintf('Outlier 2: <qU> = %.1f eV,  all std = %.0f meV , after %s std = %.0f meV \n',...
%     d.A.RunData.qU(8),1e3.*std(d.A.SingleRunData.qU(8,:)),...
%     d.A.SingleRunData.StartTimeStamp(46),1e3.*std(d.A.SingleRunData.qU(8,46:end)));





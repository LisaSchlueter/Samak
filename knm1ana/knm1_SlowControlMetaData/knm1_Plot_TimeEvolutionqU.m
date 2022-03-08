% plot mean retarding potential for knm1
% for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);
FontSize = 22;
if exist(savefile,'file')
    load(savefile);
else
    return
end

GetFigure;
s1 = subplot(4,1,[1:3]);
plot(linspace(18400,18640,10),zeros(10,1),'-','LineWidth',2,'Color',rgb('Silver'));
hold on;
for i=1:39
dscatter(1e-10.*rand(274,1)+repmat(R.RunData.qU(i),1,274)',...
    1e3*(R.SingleRunData.qU(i,:)-mean(R.SingleRunData.qU(i,:)))',...
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


qU = R.SingleRunData.qU;

s2 = subplot(4,1,4);
bar(R.RunData.qU,1e3*std(qU-mean(qU,2),0,2),'FaceColor',rgb('Orange'),'EdgeColor',rgb('Orange'));
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

%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
pltfile = sprintf('%sknm1_qU_TimeEvolution.pdf',pltdir);
export_fig(pltfile);

% LiveTime = hours(R.SingleRunData.StartTimeStamp-R.SingleRunData.StartTimeStamp(1));
% x = LiveTime;
% qUrel = qU-mean(qU,2);
% 
% plot(x,1e3*qUrel)


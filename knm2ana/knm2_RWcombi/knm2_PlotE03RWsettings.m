savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename3P = [savedir,'knm2_FitE0_3Periods.mat'];

if ~exist(savename3P,'file')
    fprintf('file %s doesnt exist \n',savename3P);
    return
end

d = importdata(savename3P); % load fit results: 1 stacked fits per RW period

figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
x = 1:numel(d.E0);

y = d.E0-d.meanE0;
yErr = d.E0Err;
p1 = bar(1,y(1),'FaceColor',rgb('IndianRed'),'FaceAlpha',0.5);
hold on
p2 = bar(2,y(2),'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.5);
p3 = bar(3,y(3),'FaceColor',rgb('GoldenRod'),'FaceAlpha',0.5);
e1 = errorbar(x(1),y(1),yErr(1),'o','LineWidth',2,'Color','k','MarkerFaceColor','k');
e1.CapSize = 0;
e1 = errorbar(x(2),y(2),yErr(2),'o','LineWidth',2,'Color','k','MarkerFaceColor','k');
e1.CapSize = 0;
e1 = errorbar(x(3),y(3),yErr(3),'o','LineWidth',2,'Color','k','MarkerFaceColor','k');
e1.CapSize = 0;
hold off
PrettyFigureFormat('FontSize',24)
xlim([0.8,3.2])
xticks(x);
xticklabels({'RW 1','RW 2','RW 3'});
ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)'));
set(gca,'XMinorTick','off');
% save
saveplot = strrep(strrep(savename3P,'results','plots'),'.mat','.pdf');
export_fig(gcf,saveplot);

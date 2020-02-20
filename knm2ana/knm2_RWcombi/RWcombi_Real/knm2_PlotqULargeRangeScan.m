% plot script for knm2_qULargeRangeScan
% plot 3 RW periods together
DataType = 'Real';
freePar = 'E0 Norm';

% load data
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
d1 = importdata([savedir,sprintf('knm2_qULargeRangeScan_%s_%s_freePar-%s.mat',DataType,'KNM2_RW1',strrep(freePar,' ',''))]);
d2 = importdata([savedir,sprintf('knm2_qULargeRangeScan_%s_%s_freePar-%s.mat',DataType,'KNM2_RW2',strrep(freePar,' ',''))]);
d3 = importdata([savedir,sprintf('knm2_qULargeRangeScan_%s_%s_freePar-%s.mat',DataType,'KNM2_RW3',strrep(freePar,' ',''))]);

meanGlobal = mean([d1.par(2,:),d2.par(2,:),d3.par(2,:)]);
y1 = (d1.par(2,:)-meanGlobal).*1e3;
yErr1 = d1.err(2,:).*1e3;
y2 = (d2.par(2,:)-meanGlobal).*1e3;
yErr2 = d2.err(2,:).*1e3;
y3 = (d3.par(2,:)-meanGlobal).*1e3;
yErr3 = d3.err(2,:).*1e3;
%% plot 1: relative endpoints
f1 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.45]);

[l,a] = boundedline(linspace(d1.RunData.qU(11)-d1.Q_i,0,10),zeros(10,1),500.*ones(10,1));
a.FaceColor = rgb('LightGray'); a.FaceAlpha = 0.7;
l.LineStyle = 'none';
hold on;
l1 = plot(linspace(min(d1.upperfitrange)-2,max(d1.upperfitrange)+2,10),zeros(10,1),'--','Color',rgb('Silver'),'LineWidth',2);
hold on
e1 = errorbar(d1.upperfitrange,y1,yErr1,'-o','CapSize',0,'LineWidth',l1.LineWidth,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
e2 = errorbar(d2.upperfitrange,y2,yErr2,'-o','CapSize',0,'LineWidth',l1.LineWidth,'Color',rgb('ForestGreen'),'MarkerFaceColor',rgb('ForestGreen'));
e3 = errorbar(d3.upperfitrange,y3,yErr3,'-o','CapSize',0,'LineWidth',l1.LineWidth,'Color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));

hold off;
PrettyFigureFormat('FontSize',24);
xlabel('Upper fit boundary below E_0 (eV)');
ylabel(sprintf('{\\itE}_0 - \\langle{\\itE}_0\\rangle (meV)'))
xlim([min(d1.upperfitrange)-2,max(d1.upperfitrange)+2])
ylim([min(y3-yErr3)*1.1,max(y2+yErr2).*2.2])
leg = legend([e1,e2,e3,a],'RW setting 1','RW setting 2','RW setting 3',sprintf(' Data used for m_\\nu fit'));%,sprintf('\\langle{\\itE}_0\\rangle'));
t = title(sprintf('Free fit parameter: %s , Lower fit boundary: -90 eV',d1.RunArg{12}));
t.FontWeight = 'normal';
t.FontSize = get(gca,'FontSize');
leg.Location='northeast';
legend boxon;
leg.EdgeColor = rgb('Silver');
leg.Color = rgb('White');

% save plot
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/plots/'];
plotname = [plotdir,sprintf('knm2_qULargeRangeScan_%s_all3RW.pdf',DataType)];
export_fig(f1,plotname);
fprintf('save plot to file %s \n',plotname);

%% plot difference with respect to RW setting 1 and with respect to largest range (-90 to -10)
f2 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.45]);

[l,a] = boundedline(linspace(d1.RunData.qU(11)-d1.Q_i,0,10),zeros(10,1),500.*ones(10,1));
a.FaceColor = rgb('LightGray'); a.FaceAlpha = 0.7;
l.LineStyle = 'none';
hold on;
e1 = errorbar(d1.upperfitrange,y1-y1-(y1(end)-y1(end)),yErr1,'-o','CapSize',0,'LineWidth',2,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
e2 = errorbar(d2.upperfitrange,y2-y1-(y2(end)-y1(end)),yErr2,'-o','CapSize',0,'LineWidth',2,'Color',rgb('ForestGreen'),'MarkerFaceColor',rgb('ForestGreen'));
e3 = errorbar(d3.upperfitrange,y3-y1-(y3(end)-y1(end)),yErr3,'-o','CapSize',0,'LineWidth',2,'Color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));

hold off;
PrettyFigureFormat('FontSize',24);
xlabel('Upper fit boundary below E_0 (eV)');
ylabel(sprintf('\\Delta{\\itE}_0 - \\Delta{\\itE}_0 ([-90,-10] eV) (meV)'))
xlim([min(d1.upperfitrange)-2,max(d1.upperfitrange)+2])
ylim([min(-yErr3)*1.2,max(y2-y1-(y2(end)-y1(end))+yErr2).*1.1])
leg = legend([e1,e2,e3,a],'RW setting 1','RW setting 2','RW setting 3',sprintf(' Data used for m_\\nu fit'));%,sprintf('\\langle{\\itE}_0\\rangle'));
t = title(sprintf('Free fit parameter: %s , Lower fit boundary: -90 eV',d1.RunArg{12}));
t.FontWeight = 'normal';
t.FontSize = get(gca,'FontSize');
leg.Location='northeast';
legend boxon;
leg.EdgeColor = rgb('Silver');
leg.Color = rgb('White');

% save plot
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/plots/'];
plotname = [plotdir,sprintf('knm2_qULargeRangeScan_%s_all3RW_diff.pdf',DataType)];
export_fig(f2,plotname);
fprintf('save plot to file %s \n',plotname);

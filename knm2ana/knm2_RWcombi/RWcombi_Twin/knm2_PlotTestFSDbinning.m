%-------------------------------------------------------------------------------------
% plot results from knm2_TestFSDbinning for 2 different TeStep together in one plot
% this script is only doing the plot, no calculation
% to get files run knm2_TestFSDbinning.m
% Lisa, December 19
%-------------------------------------------------------------------------------------
RunList = 'KNM2_Prompt';
range = 40;

savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename1 = [savedir,sprintf('knm2_TestFSDbinning_%s_%.0feVrange_%.2geV-TeBin.mat',RunList,range,0.05)];
savename2 = [savedir,sprintf('knm2_TestFSDbinning_%s_%.0feVrange_%.2geV-TeBin.mat',RunList,range,0.1)];
savename3 = [savedir,sprintf('knm2_TestFSDbinning_%s_%.0feVrange_%.2geV-TeBin.mat',RunList,range,0.15)];
savename4 = [savedir,sprintf('knm2_TestFSDbinning_%s_%.0feVrange_%.2geV-TeBin.mat',RunList,range,0.20)];

d1 = importdata(savename1,'file');
d2 = importdata(savename2,'file');
d3 = importdata(savename3,'file');
d4 = importdata(savename4,'file');
%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 = plot(d1.BinningFactorAll,d1.mNuSqShift*1e3,'o','LineWidth',2,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
hold on;
p2 = plot(d2.BinningFactorAll,d2.mNuSqShift*1e3,'o','LineWidth',2,'Color',rgb('FireBrick'),'MarkerFaceColor',rgb('FireBrick'));
p3 = plot(d3.BinningFactorAll,d3.mNuSqShift*1e3,'o','LineWidth',2,'Color',rgb('GoldenRod'),'MarkerFaceColor',rgb('GoldenRod'));
p4 = plot(d4.BinningFactorAll,d4.mNuSqShift*1e3,'o','LineWidth',2,'Color',rgb('ForestGreen'),'MarkerFaceColor',rgb('ForestGreen'));
PrettyFigureFormat('FontSize',24);
xlabel('Bin enhancement factor');
ylabel(sprintf('\\Delta m_\\nu^2 (meV^2)'));
xticks(d1.BinningFactorAll);
xlim([d1.BinningFactorAll(1)-0.2,d1.BinningFactorAll(end)+0.2]);
leg = legend([p1,p2,p3,p4],sprintf('\\Delta E = %0.0f meV',d1.TeStep*1e3),...
    sprintf('\\Delta E = %0.0f meV',d2.TeStep*1e3),sprintf('\\Delta E = %0.0f meV',d3.TeStep*1e3),...
    sprintf('\\Delta E = %0.0f meV',d4.TeStep*1e3));
legend boxoff
leg.FontSize = get(gca,'FontSize');
set(gca,'XMinorTick','off');

savenameplot = strrep(strrep(savename1,'results','plots'),'.mat','_TeSteps.pdf');
export_fig(f1,savenameplot);
fprintf('save plot to %s \n',savenameplot);



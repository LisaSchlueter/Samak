% knm2 press release flash talk
% chi2-profile: numass uncertainty

savename = [getenv('SamakPath'),'tritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_Real_UniformScan_mNu_Knm2_UniformFPD_chi2CMShape_SysBudget40_NP1.112_FitParE0BkgNorm_nFit50_KNM2_0p1eV_min-2.5_max2.5.mat'];
d = importdata(savename);

Language = 'eng';

mNuSq = reshape(d.ScanResults.ParScan,100,1);
chi2 = reshape(d.ScanResults.chi2min,100,1);
mNuSq(1) = [];
chi2(1) = [];

Idx_bf = find(chi2==min(chi2));

AxisColor  = 'k';
pltdir = [getenv('SamakPath'),'EasyVis/plots/'];
pltname = sprintf('%sProfileChi2_%s',pltdir,AxisColor);
if strcmp(Language,'eng')
    pltname = [pltname,'_eng'];
end
MakeDir(pltdir);
%% plot
chi2min = min(chi2);

%mNu = sqrt(mNuSq(mNuSq>=0));
%chi2_p = chi2(mNuSq>=0);

GetFigure
b1 = bar(mNuSq(mNuSq>0.8),chi2(mNuSq>0.8),'FaceAlpha',1,'FaceColor',rgb('Red'),'EdgeColor','none','BarWidth',0.5);
hold on;
b2 = bar(mNuSq(mNuSq<=0.8),chi2(mNuSq<=0.8),'FaceAlpha',1,'FaceColor',rgb('Red'),'EdgeColor','none','BarWidth',0.5);

%bar(mNuSq(Idx_bf),chi2(Idx_bf)-min(chi2)+1,'FaceAlpha',1,'FaceColor',rgb('DarkOrange'),'EdgeColor','none','BarWidth',0.05);
PrettyFigureFormat('FontSize',28);
if ~strcmp(AxisColor,'k')
    set(gcf,'Color','none');
    set(gca,'Color','none');
end
if strcmp(Language,'eng')
    xlabel('Neutrino mass')
    ylabel('Difference');
else
    xlabel('Neutrinomasse')
    ylabel('Unterschied');
end
yticks([]);
xticks([0 1 2]);
xticklabels({'0 eV','1 eV','2 eV'});
xlim([0 2])
ylim([26.5 54]);
set(gca,'XMinorTick','off');
ax = gca;
ax.YAxis.Color = AxisColor;
ax.XAxis.Color = AxisColor;
box off
if ~strcmp(AxisColor,'k')
    export_fig(gcf,[pltname,'.png'],'-png','-transparent','-m10');
else
    export_fig(gcf,[pltname,'.png'],'-png','-m10');
end

bg = bar(mNuSq(6),chi2(6),'FaceAlpha',1,'FaceColor',rgb('LimeGreen'),'EdgeColor','none','BarWidth',0.025);
plot(d.ScanResults.BestFit.par(1)+0.028,chi2min,'.','MarkerSize',30,'Color',AxisColor);
if ~strcmp(AxisColor,'k')
    export_fig(gcf,[pltname,'_bf.png'],'-png','-transparent','-m10');
else
    export_fig(gcf,[pltname,'_bf.png'],'-png','-m10');
end

bg2 = bar(mNuSq(mNuSq<0.8),chi2(mNuSq<0.8),'FaceAlpha',1,'FaceColor',rgb('LimeGreen'),'EdgeColor','none','BarWidth',0.5);
plot(d.ScanResults.BestFit.par(1)+0.028,chi2min,'.','MarkerSize',30,'Color',AxisColor);
if ~strcmp(AxisColor,'k')
    export_fig(gcf,[pltname,'_UpperLim.png'],'-png','-transparent','-m10');
else
    export_fig(gcf,[pltname,'_UpperLim.png'],'-png','-m10');
end
% %%
% export_fig(gcf,[pltname,'.pdf']);
% errorbar(d.ScanResults.BestFit.par(1)-0.07,34.5,0.95,'horizontal','-','Color',rgb('Black'),'LineWidth',3,'CapSize',13);
% plot(d.ScanResults.BestFit.par(1),34.5,'.','MarkerSize',30,'Color',rgb('Black'));
% export_fig(gcf,[pltname,'_err.pdf']);
% KNM2 LARA Data
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefileCorrMat = sprintf('%sknm1_LARAData.mat',savedir);
savefileCovMat = sprintf('%sknm1_LARADataCovMat_5000samples.mat',savedir);

dCov = importdata(savefileCovMat);
dCorr = importdata(savefileCorrMat);


%%
GetFigure;
corplot(dCorr.CorrMat);

xticklabels({'T_2','HT','DT'});
yticklabels({'T_2','HT','DT'});

PrettyFigureFormat('FontSize',28);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');

for i=1:3
    for j=1:3
        if dCorr.CorrMat(i,j)>0
            text(i+0.25,j+0.5,sprintf('%.3f',dCorr.CorrMat(i,j)),'FontSize',26);
        else
            text(i+0.2,j+0.5,sprintf('%.3f',dCorr.CorrMat(i,j)),'FontSize',26);
        end
    end
end

c = colorbar;
c.Label.String = sprintf('Correlation coefficient');
c.Label.FontSize = 24;
ax.XAxisLocation = 'bottom';
pltdir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/plots/'];
MakeDir(pltdir)
pltname = sprintf('knm1_LARACorrMat.pdf',pltdir);
export_fig(pltname);

% display TASR correlation
% plot for PhD thesis

% look at variations of source activity within a scan
% correlation

savedir = [getenv('SamakPath'),'knm2ana/knm2_TASR/results/'];
savefile = sprintf('%sknm2_RunObj.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

[WGTS_TASR_RelErr,SubRunActivity,TASR_CorrMat] = R.Get_DataDriven_RelErr_TASR;

BkgIdxStart = find(R.RunData.qU>R.ModelObj.Q,1);
qU = R.RunData.qU(R.exclDataStart:BkgIdxStart);
CorrMat = TASR_CorrMat(R.exclDataStart:BkgIdxStart,R.exclDataStart:BkgIdxStart);
%% 
f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
imagesc(CorrMat);
pbaspect([1 1 1]);
PrettyFigureFormat('FontSize',18);

% x,y axis
xticks([1:1:numel(qU)+0.5]);
yticks([1:1:numel(qU)+0.5]);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
EmptyStrings = strings(numel(qU)-2,1);
xticklabels({sprintf('%.0f eV',qU(1)),EmptyStrings{:},sprintf('%.0f eV',qU(end))});
yticklabels('');
xlabel('Retarding energy bin');
ylabel('Retarding energy bin');
ax = gca;
ax.YLabel.Position(1) = -0.4;
ax.XLabel.Position(2) = 25.5;

% colorbar
cb = colorbar;
cb.Label.String = sprintf('Correlation coefficient');
cb.Label.FontSize = ax.XLabel.FontSize;
cb.Position(1) = 0.86;
set(cb,'LineWidth',1.5);
colormap(flipud(hot));

%%
pltdir = [getenv('SamakPath'),'knm2ana/knm2_TASR/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm2_TASR_CorrMat.pdf',pltdir);
export_fig(pltname);
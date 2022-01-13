% display TASR systematics
% plot for PhD thesis

% look at variations of source activity within a scan
% look at variations of column density within a scan
% plot for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

[WGTS_TASR_RelErr,SubRunActivity,TASR_CorrMat] = R.Get_DataDriven_RelErr_TASR;

BkgIdxStart = find(R.RunData.qU>R.ModelObj.Q,1);

SubRunActivity = SubRunActivity(R.exclDataStart:BkgIdxStart,:);
MeanActivity = mean(mean(SubRunActivity));

qU = R.RunData.qU(R.exclDataStart:BkgIdxStart,:);
qU_singlerun = R.SingleRunData.qU(R.exclDataStart:BkgIdxStart,:);

SubRunActivity_scatter = reshape((SubRunActivity-MeanActivity)./MeanActivity,1,R.nRuns.*numel(qU))';
qU_scatter = repmat(qU,R.nRuns,1);

%% test 

%GetFigure
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
% relative standard deviation
[pl, ps] = boundedline(qU,(mean(SubRunActivity,2)-MeanActivity)./MeanActivity,std(SubRunActivity,0,2)./MeanActivity);
ps.FaceColor = rgb('LightSalmon'); ps.EdgeColor = rgb('LightSalmon'); pl.LineStyle = 'none';
hold on;

% relative error of the mean
[pl2, peom] = boundedline(qU,(mean(SubRunActivity,2)-MeanActivity)./MeanActivity,std(SubRunActivity,0,2)./(MeanActivity.*sqrt(R.nRuns)));
peom.FaceColor = rgb('Red'); peom.EdgeColor = rgb('Red'); pl2.LineStyle = 'none';

for i=1:numel(qU)
    % activity distribution for each scan-step
    % relative to mean activit of this scan-step
    sp = dscatter(qU_singlerun(i,:)',(SubRunActivity(i,:)'-MeanActivity)./MeanActivity,'Msize',20); % density scatter plot
    hold on;
end
hold on;

%plot(qU,zeros(size(mean(SubRunActivity,2))),'-r','LineWidth',2);%zeros(size(mean(SubRunActivity,2)))
xlabel('Retarding energy');
ylabel('Relative source activity');%sprintf('[\\rhod - \\langle\\rhod\\rangle]/\\langle\\rhod\\rangle (%%)'));
PrettyFigureFormat('FontSize',18);
c = colorbar;
c.Label.String = 'Scan density'; 
c.Label.FontSize = get(gca,'FontSize')+2;
colormap(flipud(colormap('winter')));
ax = gca;
ax.XAxis.Exponent = 0;
xlim([min(qU)-2 max(qU)+2]);
xticks(18535:10:18575);

leg = legend([ps,peom],'Standard deviation','Error of the mean','Location','southeast');
legend box off;

pltdir = [getenv('SamakPath'),'knm1ana/knm1_TASR/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm1_TASR.pdf',pltdir);
export_fig(pltname);

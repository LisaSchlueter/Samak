% display TASR systematics
% plot for PhD thesis
% knm2 standard settings for uniform FPD (stat. only)

% look at variations of source activity within a scan
% look at variations of column density within a scan
% plot for PhD thesis
savedir = [getenv('SamakPath'),'knm2ana/knm2_TASR/results/'];
savefile = sprintf('%sknm2_RunObj.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else     
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Real',...
        'fixPar','mNu E0 Norm Bkg',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',40,...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'RingMerge','None',...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    R = MultiRunAnalysis(RunAnaArg{:});
    %%
    R.exclDataStart = R.GetexclDataStart(40);
    R.InitModelObj_Norm_BKG;
    R.Fit;
    MakeDir(savedir);
    save(savefile,'R');
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

GetFigure

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


%plot(qU,zeros(size(mean(SubRunActivity,2))),'-r','LineWidth',2);%zeros(size(mean(SubRunActivity,2)))
xlabel('Retarding energy');
ylabel('Relative source activity');%sprintf('[\\rhod - \\langle\\rhod\\rangle]/\\langle\\rhod\\rangle (%%)'));
PrettyFigureFormat('FontSize',22);
c = colorbar;
c.Label.String = 'Scan density'; 
c.Label.FontSize = get(gca,'FontSize')+2;
colormap(flipud(colormap('winter')));
ax = gca;
ax.XAxis.Exponent = 0;
xlim([min(qU)-2 max(qU)+2]);
xticks(18534:10:18574);

leg = legend([ps,peom],'Standard deviation','Error of the mean','Location','southwest');
legend box off;
%%
pltdir = [getenv('SamakPath'),'knm2ana/knm2_TASR/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm2_TASR.pdf',pltdir);
export_fig(pltname);

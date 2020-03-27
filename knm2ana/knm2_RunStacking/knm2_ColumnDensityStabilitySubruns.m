RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'ROIFlag','14keV',...
    'SysBudget',33};

A = MultiRunAnalysis(RunAnaArg{:});
%%
qU = A.RunData.qU-18574;
WGTS_subrun     = A.SingleRunData.WGTS_CD_MolPerCm2_SubRun; % nqU x nRun
MeanWGTS_subrun = mean(WGTS_subrun);% average for every run % nRuns
RelFluct = (WGTS_subrun-MeanWGTS_subrun)./MeanWGTS_subrun;  % nqU x nRuns
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
% plot(qU,RelFluct);
% hold on;
[l,a] = boundedline(qU,100.*mean(RelFluct,2),100.*std(RelFluct,0,2));
l.LineWidth = 2; l.Color = rgb('DodgerBlue');
a.FaceColor = rgb('PowderBlue'); a.FaceAlpha=1;
leg = legend([l,a],'Average',sprintf('1\\sigma band'),'EdgeColor',rgb('Silver'),'Location','northwest');
xlim([-90 5]);
PrettyFigureFormat('FontSize',24);
xlabel('Retarding potential - 18574 (eV)');
ylabel('Rel. column density diff. (%)');
title(sprintf('KNM2 %.0f golden runs',A.nRuns),'FontWeight','normal','FontSize',get(gca,'FontSize'));
savename = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/plots/knm2_ColumnDensityStabilitySubruns.pdf'];
export_fig(f2,savename);
%%
fprintf('Average std from subrun to subrun = %.2f %% \n',mean(std(WGTS_subrun)./MeanWGTS_subrun)*100);
fprintf('Average std from run to run = %.2f %% \n',100*std(A.SingleRunData.WGTS_CD_MolPerCm2)./mean(A.SingleRunData.WGTS_CD_MolPerCm2))
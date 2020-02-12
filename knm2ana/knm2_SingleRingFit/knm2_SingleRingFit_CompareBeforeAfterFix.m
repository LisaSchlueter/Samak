% uniform fit on knm2 stacked data
% settings
RunList = 'KNM2_afterFix';
RunList2 = 'KNM2_beforeFix';
fixPar = 'E0 Norm Bkg'; % free parameter
DataType = 'Real';
FSDFlag = 'BlindingKNM2';
NonPoissonScaleFactor = 1.0;%12;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 1; % 11 == 40 eV range
chi2 = 'chi2Stat';
RunAnaArg = {'RunList',RunList,'fixPar',fixPar,'DataType',DataType,...
            'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,'exclDataStart',exclDataStart,...
            'AnaFlag',AnaFlag,'chi2',chi2,...
            'RingMerge','Full',...
            'fitter','minuit'};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});
R = RingAnalysis('RunAnaObj',A,'RingList',1:4);
R.FitRings;

%% RunAnalysis object for all runs before the fix
RunAnaArg2 = {'RunList',RunList2,'fixPar',fixPar,'DataType',DataType,...
    'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,'exclDataStart',exclDataStart,...
    'AnaFlag',AnaFlag,'chi2',chi2,...
    'RingMerge','Full',...
    'fitter','minuit'};
A2 = MultiRunAnalysis(RunAnaArg2{:});
R2 = RingAnalysis('RunAnaObj',A2,'RingList',1:4);
R2.FitRings;
%% plot
fig2 = figure('Renderer','painters');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.6,0.5]);
plot(0:5,zeros(6,1),'--','Color',rgb('Silver'),'LineWidth',2); hold on
e2= errorbar(1:4,R2.FitResult.par(:,2),R2.FitResult.err(:,2),'o',...
        'LineWidth',2,'LineStyle','none','MarkerSize',8);
    hold on;
e1 = errorbar(1:4,R.FitResult.par(:,2),R.FitResult.err(:,2),'o',...
    'LineWidth',2,'LineStyle','none','MarkerSize',8);        
xlabel('Rings')
xticks([1 2 3 4]);
xticklabels({'1,2,3','4,5,6','7,8,9','10,11,12'})
ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)'));
PrettyFigureFormat('FontSize',24);
xlim([0.5,4.5])

[linFitpar, linFiterr, linFitchi2min,linFitdof] =...
    linFit(R.RingList',R.FitResult.par(:,2),R.FitResult.err(:,2));
[linFitpar2, linFiterr2, linFitchi2min2,linFitdof2] =...
    linFit(R.RingList',R2.FitResult.par(:,2),R2.FitResult.err(:,2));

l = plot(R.RingList,(linFitpar(1).*R.RingList+linFitpar(2))','-','Color',e1.Color,'LineWidth',3);
linFitleg =  sprintf('After Fix: linear fit slope: (%.1f \\pm %.1f) meV/ring @ \\chi2 = %.1f / %.0f dof',linFitpar(1)*1e3,linFiterr(1)*1e3,linFitchi2min,linFitdof);
l2 = plot(R.RingList,(linFitpar2(1).*R.RingList+linFitpar2(2))','-','Color',e2.Color,'LineWidth',3);
linFitleg2 =  sprintf('Before Fix: linear fit slope: (%.1f \\pm %.1f) meV/ring @ \\chi2 = %.1f / %.0f dof',linFitpar2(1)*1e3,linFiterr2(1)*1e3,linFitchi2min2,linFitdof2);

leg = legend([l2,l],linFitleg2,linFitleg);
legend boxoff

ylim([-0.1,0.1]);
savedir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/RingFit/',A.DataSet)];
MakeDir(savedir);
savename = sprintf('%s_ringwise_%s_%.0frange_%s_Merge_%s.pdf',...
    'E0','BeforeAndAfterFix',90,A.RingMerge,A.DataType);
print([savedir,savename],'-dpng','-r350');
export_fig(gcf,[savedir,savename],'-painters');
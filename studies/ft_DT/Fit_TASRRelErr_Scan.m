
WGTS_TASR_RelErr = [1:10].*1e-3;
chi2min = zeros(numel(WGTS_TASR_RelErr),1);
pvalue  = zeros(numel(WGTS_TASR_RelErr),1);

for i=1:numel(WGTS_TASR_RelErr)
MRA = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','fixPar','1 5 6','exclDataStart',9);    
MRA.ComputeCM('WGTS_TASR_RelErr',WGTS_TASR_RelErr(i));
MRA.chi2 = 'chi2CM';
MRA.Fit;
chi2min(i) = MRA.FitResult.chi2min;
pvalue(i)  = 1-chi2cdf(chi2min(i),MRA.FitResult.dof);
end
%%
s1 = subplot(2,1,1);
plot(WGTS_TASR_RelErr*1e2,chi2min,'LineWidth',3,'Color',rgb('CadetBlue'));
hold on;
plot(WGTS_TASR_RelErr*1e2,MRA.FitResult.dof.*ones(numel(WGTS_TASR_RelErr),1),'--k','LineWidth',3,'Color',rgb('Gray'));
xlabel('Tritium Subrun to Subrun Activity Fluctuation (%)');
ylabel(['\chi^2 (',num2str(MRA.FitResult.dof),' dof)']);
title(sprintf('Samak Fit to Stacked Runs (%s) with Systematics - %.0f eV below E0',MRA.StackFileName, MRA.ModelObj.Q_i-MRA.ModelObj.qU(9)));
PrettyFigureFormat;
set(gca,'FontSize',18);
hold off;

s2 = subplot(2,1,2);
plot(WGTS_TASR_RelErr*1e2,pvalue,'LineWidth',3,'Color',rgb('CadetBlue'));
hold on;
plot(WGTS_TASR_RelErr*1e2,0.05.*ones(numel(WGTS_TASR_RelErr),1),'--k','LineWidth',3,'Color',rgb('Gray'));
xlabel('Tritium Subrun to Subrun Activity Fluctuation (%)');
ylabel('p-value');
PrettyFigureFormat;
set(gca,'FontSize',18);
set(gca,'YScale','log');
hold off;
xlim([min(WGTS_TASR_RelErr*1e2) max(WGTS_TASR_RelErr*1e2)]);

linkaxes([s1,s2], 'x');
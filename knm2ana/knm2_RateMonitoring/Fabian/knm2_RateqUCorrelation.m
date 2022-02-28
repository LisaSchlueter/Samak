nSamples = 1e4;
nRuns = 361;
TBDIS_i = 6.6788e+04; % mean rate at rate monitor point
%%
TBDIS = TBDIS_i+sqrt(TBDIS_i).*randn(nRuns,nSamples);
time = ((1:nRuns).*ones(nSamples,1))';

%% correlation coefficients
CorrCoeff = corr((1:nRuns)',TBDIS);
GetFigure;
h1 = histogram(CorrCoeff,'FaceColor',rgb('PowderBlue'),'FaceAlpha',1);
xlabel('Correlation coefficient (time , -300 eV rate)')
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)
savepath =  [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
MakeDir(savepath);
savename = sprintf('%sknm2_RandomRateTimeCorrelation.png',savepath);
print(savename,'-dpng','-r500');


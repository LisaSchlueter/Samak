% Create MultiRunAnalysis Object
if ~exist('SPR_R','var') || ~exist('SPR_T','var')
    [SPR_R , ~] = Knm1RealTwin_Create();
end

% Data or Twin
SPRG = SPR_R;

SPRG.FitRunList('Recompute','ON')

%% Get Single Run Model Objects
if isempty(SPRG.SingleRunObj)
SPRG.LoadSingleRunObj;
end

%% Extract Cum Frac
SR_CumFrac      = cellfun(@(x) x.CumFrac,SPRG.SingleRunObj,'UniformOutput', true);
NormFactorTBDDS = cellfun(@(x) x.NormFactorTBDDS,SPRG.SingleRunObj,'UniformOutput', true);
Tpurity         = cellfun(@(x) x.WGTS_epsT,SPRG.SingleRunObj,'UniformOutput', true);


%% NormFit Versus Other Fit Parameters
NormFitSR = figure(1);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
subplot(2,1,1)
scatter(SPRG.SingleRun_FitResults.chi2Stat.E0,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
xlabel('Single Run Fit E_0 (eV)');
ylabel('Single Run Fit Normalization (a.u.)');
PrettyFigureFormat
subplot(2,1,2)
scatter(SPRG.SingleRun_FitResults.chi2Stat.B,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
xlabel('Single Run Fit Background (cps)');
ylabel('Single Run Fit Normalization (eV)');
PrettyFigureFormat
export_fig(1,'plots/knm1_singlerun_normfit_1.png','-r300');

%% Plot  NormFit VS qUmin
NormFitSR = figure(2);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(min(SPRG.SingleRunData.qU),SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Run-wise fit normalization');
xlabel('min(qU)');
PrettyFigureFormat
export_fig(2,'plots/knm1_singlerun_normfit_2.png','-r300');

%% BiPlot: Endpoint Fit VS qUmin VS N
x = SPRG.SingleRun_FitResults.chi2Stat.E0';
y = min(SPRG.SingleRunData.qU)';
z = SPRG.SingleRun_FitResults.chi2Stat.N';

xlin = linspace(min(x),max(x),50);
ylin = linspace(min(y),max(y),50);

[X,Y] = meshgrid(xlin,ylin);

f = scatteredInterpolant(x,y,z);
Z = f(X,Y);

NormFitSR = figure(3);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
imagesc(xlin-SPRG.ModelObj.Q_i,ylin-SPRG.ModelObj.Q_i,Z) %interpolated
ylabel(sprintf('qUmin - %.3f (eV)',SPRG.ModelObj.Q_i));
xlabel(sprintf('Single Run Fit E_0 - %.3f (eV)',SPRG.ModelObj.Q_i));
title('Single run Fit Normalization')
colorbar
PrettyFigureFormat
export_fig(3,'plots/knm1_singlerun_normfit_3.png','-r300');


%% Plot  NormFit Versus Run Number
NormFitSR = figure(4);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(SPRG.RunList,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Single Run Fit Normalization');
xlabel('Run Number');
PrettyFigureFormat
export_fig(4,'plots/knm1_singlerun_normfit_4.png','-r300');

%% Plot  NormFit Versus RhoD
NormFitSR = figure(5);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRun_FitResults.chi2Stat.N+1,'s','LineWidth',2);

P = polyfit(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRun_FitResults.chi2Stat.N+1,1);
yfit = P(1)*SPRG.SingleRunData.WGTS_CD_MolPerCm2+P(2);
hold on;
plot(SPRG.SingleRunData.WGTS_CD_MolPerCm2,yfit,'r-.','LineWidth',5,'Color',rgb('SteelBlue'));
hold off

ylabel('Single Run Fit Normalization');
xlabel('Column Density (mol/cm^2)');
PrettyFigureFormat
export_fig(5,'plots/knm1_singlerun_normfit_5.png','-r300');

%% Plot  NormFit VS [TT] / [HT] / [DT]
CumFracFig = figure(6);
set(CumFracFig,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=subplot(3,1,1)
scatter(SPRG.SingleRunData.WGTS_MolFrac_TT,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Fit Norm');
xlabel('[TT]');
PrettyFigureFormat
s2=subplot(3,1,2)
scatter(SPRG.SingleRunData.WGTS_MolFrac_HT,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Fit Norm');
xlabel('[HT]');
PrettyFigureFormat
s2=subplot(3,1,3)
scatter(SPRG.SingleRunData.WGTS_MolFrac_DT,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Fit Norm');
xlabel('[DT]');
PrettyFigureFormat
export_fig(6,'plots/knm1_singlerun_normfit_6.png','-r300');

%% Plot  NormFit Versus Endpoint
NormFitSR = figure(7);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(SPRG.SingleRun_FitResults.chi2Stat.E0-SPRG.ModelObj.Q_i,SPRG.SingleRun_FitResults.chi2Stat.N,'s','LineWidth',2);
ylabel('Single Run Fit Normalization');
xlabel(sprintf('Single Run Fit E_0 - %.3f (eV)',SPRG.ModelObj.Q_i));
PrettyFigureFormat
export_fig(7,'plots/knm1_singlerun_normfit_7.png','-r300');

%% Plot  RhoD Versus Endpoint
NormFitSR = figure(8);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRun_FitResults.chi2Stat.E0-SPRG.ModelObj.Q_i,'s','LineWidth',2);
xlabel('Column Density (mol/cm^2)');ylabel(sprintf('Single Run Fit E_0 - %.3f (eV)',SPRG.ModelObj.Q_i));
PrettyFigureFormat
export_fig(8,'plots/knm1_singlerun_normfit_8.png','-r300');

%% BiPlot: Endpoint Fit VS rhoD VS N
x = SPRG.SingleRun_FitResults.chi2Stat.E0';
y = SPRG.SingleRunData.WGTS_CD_MolPerCm2';
z = SPRG.SingleRun_FitResults.chi2Stat.N';

xlin = linspace(min(x),max(x),25);
ylin = linspace(min(y),max(y),25);

[X,Y] = meshgrid(xlin,ylin);

f = scatteredInterpolant(x,y,z);
Z = f(X,Y);

NormFitSR = figure(9);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
imagesc(xlin-SPRG.ModelObj.Q_i,ylin,Z) %interpolated
ylabel('\rho d (mol/cm^2)');
xlabel(sprintf('Single Run Fit E_0 - %.3f (eV)',SPRG.ModelObj.Q_i));
title('Single run Fit Normalization')
h = colorbar;
set(get(h,'title'),'string','N-1');
PrettyFigureFormat
export_fig(9,'plots/knm1_singlerun_normfit_9.png','-r300');

%%
%% Plot  RhoD Versus Endpoint
P(1)=8.273976836242512e-18;
P(2)=0.141530047843744;
NormFitSR = figure(10);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,NormFactorTBDDS');
hold on
s2=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,NormFactorTBDDS'.*(P(1)*SPRG.SingleRunData.WGTS_CD_MolPerCm2+P(2)));
hold off
legend([s1 s2],'Not Corrected',' Corrected');
PrettyFigureFormat
export_fig(10,'plots/knm1_singlerun_normfit_10.png','-r300');


%% Plot  RhoD Versus Endpoint
P(1)=8.273976836242512e-18;
P(2)=0.141530047843744;
NormFitSR = figure(11);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRunData.WGTS_CD_MolPerCm2);
hold on
s2=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRunData.WGTS_CD_MolPerCm2.*(P(1)*SPRG.SingleRunData.WGTS_CD_MolPerCm2+P(2)));
hold off
legend([s1 s2],'Not Corrected',' Corrected');
PrettyFigureFormat

%% Plot  RhoD Versus Endpoint
P(1)=8.273976836242512e-18;
P(2)=0.141530047843744;
NormFitSR = figure(11);
set(NormFitSR,'Units', 'normalized','Position', [0.1, 0.1, 0.8, 0.8]);
s1=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,SPRG.SingleRunData.WGTS_CD_MolPerCm2./SPRG.SingleRunData.WGTS_CD_MolPerCm2);
hold on
s2=scatter(SPRG.SingleRunData.WGTS_CD_MolPerCm2,(P(1)*SPRG.SingleRunData.WGTS_CD_MolPerCm2+P(2)));
hold off
legend([s1 s2],'Not Corrected',' Corrected');
PrettyFigureFormat
export_fig(11,'plots/knm1_singlerun_normfit_11.png','-r300');




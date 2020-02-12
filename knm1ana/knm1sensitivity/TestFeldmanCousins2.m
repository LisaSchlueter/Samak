% Feld man Counsins Belt 
% MC method: fit statistically fluctuated spectra
% calculate Delta chi2 form fit
% construct confidence interval with feldman cousins ordering principle

%% set up model
D = MultiRunAnalysis('RunList','KNM1','FSDFlag','Sibille0p5eV','fitter',...
    'minuit','minuitOpt','min;migrad','exclDataStart',14,'chi2','chi2Stat');
S = RunSensitivity('RunAnaObj',D);
%% calculate / get x1,x2
 S.ComputeFC('mNuSq_t',[0,0.1,0.15,0.2:0.1:0.7,0.9,1.2,1.35,1.5],'nSamples',3000);
 S.FC_ComputeX1X2;
%% plot
f111 = figure('Renderer','opengl');
set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
PlotArg = {'-o','Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',3,'MarkerSize',8};

a1 = area([S.FC_x1;linspace(S.FC_x1(end),3,1)],[S.FC_mNuSqTrue;2],'FaceColor',rgb('SteelBlue'),'FaceAlpha',0.1,'LineStyle','none');
hold on;
a2 =area(S.FC_x2,S.FC_mNuSqTrue,0,'FaceColor',rgb('White'),'FaceAlpha',1,'LineStyle','none');
p1 =plot(S.FC_x1,S.FC_mNuSqTrue,PlotArg{:});
p2 =plot(S.FC_x2,S.FC_mNuSqTrue,PlotArg{:});
PrettyFigureFormat;
xlabel(sprintf('measured m_\\beta^2 (eV^2)'));
ylabel(sprintf('true m_\\beta^2 (eV^2)'));
xlim([min(S.FC_x1), max(S.FC_x2)]);
ylim([0,max(S.FC_mNuSqTrue)]);
mNuMeasured = -0.98;
x1 = S.FC_x1(~isnan(S.FC_x1));
yticks(0:0.2:max(S.FC_mNuSqTrue))
mNuLimit = interp1(x1,S.FC_mNuSqTrue(~isnan(S.FC_x1)),mNuMeasured,'spline');
plot(mNuMeasured*[1,1],[0,mNuLimit],'--','LineWidth',3.5,'color',rgb('Orange'));
plot([min(S.FC_x1),mNuMeasured],[mNuLimit,mNuLimit],'--','LineWidth',3.5,'color',rgb('Orange'));
leg = legend(p2,'90% C.L. (stat only)','Location','northwest');
legend boxoff;
%%
savedir = [getenv('SamakPath'),'/knm1ana/knm1_unblinding/plots/'];
savefile = [savedir,'FCbelt.png'];
print(f111,savefile,'-dpng','-r450');



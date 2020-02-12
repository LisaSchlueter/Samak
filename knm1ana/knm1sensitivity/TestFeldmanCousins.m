M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11');
S = RunSensitivity('RunAnaObj',M);
S.ComputeFCR('nSamples',50,'mNuSq_t',0.5);%[-1.5:0.5:-0.5,0,0.1,0.2,0.3:0.2:0.9,2,5:0.5:6]);
%% plot delta chi2 
 f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

mNuTest = 1;
x = squeeze(S.NFitResults.par(mNuTest,1,:));
DeltaChi2 = S.FC_Chi2TrueM(mNuTest,:)-S.FC_Chi2Mbest(mNuTest,:);
plot(x,DeltaChi2,'x','Color',rgb('GoldenRod'));
hold on;
plot(linspace(min(x),max(x),10),prctile(DeltaChi2,90).*ones(10,1),...
    '--','Color',rgb('FireBrick'),'LineWidth',3);
hold off;
PrettyFigureFormat;
%xlim([-2,2]);
%ylim([-0.1,6
xlabel(sprintf('m^2_\\nu (eV^2)'));
ylabel(sprintf('\\Delta \\chi^2'));
leg = legend('toy experiments',sprintf('\\chi^2_c')); legend boxoff;
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/plots/'];
leg.Location = 'northwest';
savename = [savedir,'FC_construction.png'];
set(gca,'FontSize',24);
%print(savename,'-dpng','-r450')
%% plot neutrino masses
 f2 = figure('Renderer','opengl');
set(f2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);
mNuTest = 1;
h1 = histogram(x);
xlabel(sprintf('m^2_\\nu (eV^2)'));
ylabel('toy experiments');
PrettyFigureFormat;
h1.FaceColor = rgb('GoldenRod'); h1.FaceAlpha = 0.9;
leg = legend(sprintf('m^2_\\nu (model) = %.1f eV^2',S.FC_mNuSq(mNuTest)));
leg.Location = 'northwest';
legend boxoff;
set(gca,'FontSize',24);
savename = [savedir,'FC_mNuDist.png'];
print(f2,savename,'-dpng','-r450')
%hold on;
%% cl beld
 f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

%low = S.FC_LimitLow;
%low(low==0) = NaN;
plot(S.FC_mNuSq,S.FC_LimitLow,'-o','LineWidth',3,'Color',rgb('FireBrick'));
hold on;
plot(S.FC_mNuSq,S.FC_LimitUp,'-o','LineWidth',3,'Color',rgb('FireBrick'));
hold off;
xlabel(sprintf('m^2_\\nu measured (eV^2)'));
ylabel(sprintf('F.C. confidence interval %.0f C.L. (eV^2)',S.ConfLevel*100));
PrettyFigureFormat;
savename = [savedir,'FC_preliminary.png'];
print(f1,savename,'-dpng','-r450')
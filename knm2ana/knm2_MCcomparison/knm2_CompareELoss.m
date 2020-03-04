
KaFit = 'OFF'; % with or without Kafit
RecomputeFlag = 'ON';
%% Load samak
maxE = 9288;
ELossBinStep = 0.01;
minE=-maxE; NbinE = (maxE-minE)/ELossBinStep;
E = minE:ELossBinStep:maxE;
EIndex = (E>=0 & E<90);

savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
savename = sprintf('%sELoss_KatrinT2_%.2feVbinning.mat',savedir,ELossBinStep);
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    A = ref_FakeRun_KNM2_RFcomparison;
    A.recomputeRF = 'ON';
    [~, ElossFunctions] = A.ComputeELossFunction('E',E);
     save(savename,'ElossFunctions','E','ELossBinStep');
end

Es = E(EIndex);
S1 = ElossFunctions(1,EIndex);            % one scattering
S2 = ElossFunctions(2,EIndex).*ELossBinStep;     % two scatterings
S3 = ElossFunctions(3,EIndex).*ELossBinStep.^2;  % three scatterings

%% load Fitrium
savenameF = sprintf('%sELoss_KatrinT2_Fitrium.dat',savedir);
d = importdata(savenameF);
Ef = d(:,1)';
F1 = d(:,2)';
F2 = d(:,3)';
F3 = d(:,4)';

%% load kafit
savenameK = sprintf('%sELoss_KatrinT2_KaFit.dat',savedir);
dk = importdata(savenameK);
Ek = dk(:,1)';
StartIndex = find(Ek==0);
Ek = dk(StartIndex:end-1,1)';
K1 = dk(StartIndex:end-1,2)';

%% difference
fig1 = figure('Renderer','opengl');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6 ,0.6]);
switch KaFit
    case 'ON'
        p1 = plot(Es,S1-F1,'-.','Linewidth',3,'Color',rgb('GoldenRod'));
        hold on;
        p2 = plot(Es,S1-K1,'-','Linewidth',3,'Color',rgb('DodgerBlue'));
        p3 = plot(Es,F1-K1,':','Linewidth',3,'Color',rgb('IndianRed'));
        leg = legend([p2,p1,p3],'Samak - KaFit','Samak - Fitrium','Fitrium - KaFit');
        leg.Title.String = 'First scattering';
    case 'OFF'
        p1 = plot(Es,S1-F1,'-','Linewidth',3,'Color',rgb('FireBrick'));
        hold on;
        p2 = plot(Es,S2-F2,'Linewidth',3,'Color',rgb('GoldenRod'),'LineStyle',':');
        p3 = plot(Es,S3-F3,'Linewidth',3,'Color',rgb('RoyalBlue'),'LineStyle','-.');
        leg = legend('1 scattering', '2 scatterings','3 scatterings');
        leg.Title.String = 'Samak - Fitrium';
        leg.Location = 'southeast';
        hold off;
end
PrettyFigureFormat;
set(gca,'FontSize',24);
xlim([0 90]);
xlabel(['energy loss ',char(949),' (eV)']);
ylabel(['probability f(',char(949),') diff.']);
leg.EdgeColor = rgb('Silver');
leg.FontSize = 22;

plotdir = strrep(savedir,'results','plots');
switch KaFit
    case 'ON'
        savename = sprintf('%sELoss_Diff',plotdir);
    case 'OFF'
        savename = sprintf('%sELoss_Diff_Fitrium',plotdir);
end
export_fig(fig1,[savename,'.pdf']);
print(fig1,[savename,'.png'],'-dpng','-r300');

 
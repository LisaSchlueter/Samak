

%% Load samak
maxE = 9288;
ELossBinStep = 0.1;
minE=-maxE; NbinE = (maxE-minE)/ELossBinStep;
E = minE:ELossBinStep:maxE;
EIndex = (E>=0 & E<90);

savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
savename = sprintf('%sELoss_KatrinT2_%.2feVbinning.mat',savedir,ELossBinStep);
if exist(savename,'file') && 2==1
    load(savename);
else
    A = ref_FakeRun_KNM2_RFcomparison;
    A.recomputeRF = 'ON';
    [~, ElossFunctions] = A.ComputeELossFunction('E',E);
     save(savename,'ElossFunctions','E','Estep');
end
Es = E(EIndex);
S1 = ElossFunctions(1,EIndex);            % one scattering
S2 = ElossFunctions(2,EIndex).*Estep;     % two scatterings
S3 = ElossFunctions(3,EIndex).*Estep.^2;  % three scatterings

%% load Fitrium
savenameF = sprintf('%sELoss_KatrinT2_Fitrium.dat',savedir);
d = importdata(savenameF);
Ef = d(:,1)';
F1 = d(:,2)';
F2 = d(:,3)';
F3 = d(:,4)';

%% difference
fig5 = figure('Renderer','opengl');
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6 ,0.6]);
plot(Es,S1-F1,'Linewidth',3,'Color',rgb('FireBrick'));
hold on;
plot(Es,S2-F2,'Linewidth',3,'Color',rgb('Orange'),'LineStyle',':');
hold on;
plot(Es,S3-F3,'Linewidth',3,'Color',rgb('RoyalBlue'),'LineStyle','--');
PrettyFigureFormat;
set(gca,'FontSize',24);
xlim([0 90]);
xlabel(['energy loss ',char(949),' (eV)']);
ylabel(['probability f(',char(949),') diff.']);
leg = legend(string(1:3));
leg.EdgeColor = rgb('Silver');
leg.Location = 'southeast';
leg.FontSize = 22;
leg.Title.String = 'number of scatterings';
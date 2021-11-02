% some key parameters for knm1
% signal to background
% for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

if isempty(R.FitResult)
    R.Fit;
    save(savefile,'R','-append');
end

%% number of electrons: analysis interval
BkgRate = R.ModelObj.BKG_RateSec;
Nbkg_qU  = R.RunData.qUfrac(R.exclDataStart:end).*R.RunData.TimeSec.*BkgRate;
Nbkg     = sum(Nbkg_qU);
Nall    = sum(R.RunData.TBDIS(R.exclDataStart:end));
Nsig = Nall-Nbkg;

fprintf('---------- 40 eV interval ----------------\n');
fprintf('Total      number electrons %.2e \n',Nall);
fprintf('Signal     number electrons %.2e \n',Nsig);
fprintf('Background number electrons %.2e \n',Nbkg);

%% signal to background
Nsig_qU = R.RunData.TBDIS(R.exclDataStart:end)-Nbkg_qU;
SB_qU = Nsig_qU./Nbkg_qU;
fprintf('--------------------------------\n');
fprintf('Tot signal to background   %.2f \n',Nsig/Nbkg);
fprintf('Mean signal to background   %.2f \n',mean(SB_qU));
fprintf('Min signal to background   %.2f \n',min(SB_qU));
fprintf('Max signal to background   %.2f \n',max(SB_qU));

%%
qU = R.RunData.qU(R.exclDataStart:end);

f1 = figure('Units','normalized','Position',[1.1,1.1,0.5,0.4]);%  GetFigure;
bar(qU-R.ModelObj.Q,SB_qU,'FaceColor',rgb('Silver'),'EdgeColor',rgb('Silver'));
hold on;
p1 = plot(linspace(-40,0,27),ones(27,1),':k','LineWidth',2.5);
xlim([-40 0]);
ylim([0 170]);
xlabel(sprintf('Retarding energy - {\\itE}_0 (eV)'));
ylabel('Signal-to-background ratio');
set(gca,'YScale','log');
interp1(SB_qU(14:18),qU(14:18)-R.ModelObj.Q_i,1,'lin')
Sig2Bkg1 = interp1(SB_qU(14:18),qU(14:18)-R.ModelObj.Q_i,1,'spline');
leg = legend(p1,sprintf('Signal-to-background equilibrium at {\\itqU - E}_0 = %.0f eV',Sig2Bkg1),'Location','northeast');
legend boxoff
PrettyFigureFormat('FontSize',22);
leg.FontSize = get(gca,'FontSize');

ax = gca;
ax.XAxis.Exponent = 0;
%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
pltfile = sprintf('%sknm1_Sig2Bkg.pdf',pltdir);
export_fig(pltfile);



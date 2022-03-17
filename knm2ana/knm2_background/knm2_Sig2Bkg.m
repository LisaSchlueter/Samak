% some key parameters for knm2
% signal to background
% for PhD thesis
savedirf = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savenamef = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedirf,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');

if exist(savenamef,'file')
    d = importdata(savenamef);
    fprintf('load %s \n',savenamef);
else
    fprintf(2,'file not found %s \n',savenamef)
    return
end

R = d.A;

%% number of electrons: analysis interval based on DATA + fitted steady-state background rate 
R_i_bkg   = R.ModelObj.BKG_RateSec;
Nbkg_qU   = R.RunData.qUfrac(R.exclDataStart:end).*R.RunData.TimeSec.*R_i_bkg;
Nbkg      = sum(Nbkg_qU);
Nbkg_qU_sig  = R.RunData.qUfrac(R.exclDataStart:end-5).*R.RunData.TimeSec.*R_i_bkg; % only signal points
Nbkg_sig  = sum(Nbkg_qU_sig); % only signal points
Nall      = sum(R.RunData.TBDIS(R.exclDataStart:end));
Nsig      = Nall-Nbkg;

fprintf('---------- 40 eV interval ----------------\n');
fprintf('Total      number electrons %.2e \n',Nall);
fprintf('Signal     number electrons %.2e \n',Nsig);
fprintf('Background number electrons %.2fe6 (all points) \n',Nbkg.*1e-06);
fprintf('Background number electrons %.2fe6 (only signal points) \n',Nbkg_sig.*1e-06);
fprintf('------------------------------------------------------------\n');

% signal to background (all data points in 40 eV)
Nsig_qU = R.RunData.TBDIS(R.exclDataStart:end)-Nbkg_qU;
SB_qU = Nsig_qU./Nbkg_qU;

fprintf('--------------------------------\n');
fprintf('Tot (analysis interval) signal to background   %.2f \n',Nsig/Nbkg);
fprintf('Mean signal to background   %.2f \n',mean(SB_qU));
fprintf('Min signal to background   %.2f \n',min(SB_qU));
fprintf('Max signal to background   %.2f \n',max(SB_qU));

% signal-to-background equilibrium
qU = R.RunData.qU(R.exclDataStart:end);
Sig2Bkg1_lin = interp1(SB_qU(14:20),qU(14:20)-R.ModelObj.Q_i,1,'lin'); %cross check with linear interp (more robust, less accurate)
Sig2Bkg1 = interp1(SB_qU(14:18),qU(14:18)-R.ModelObj.Q_i,1,'spline');
fprintf('Signal to background  equilibrium  at qU-E0 = %.2f eV\n',Sig2Bkg1);


%% signal to background (only signal points in 40 eV)
Nsig_qU_sig = R.RunData.TBDIS(R.exclDataStart:end-5)-Nbkg_qU_sig;
SB_qU_sig = Nsig_qU_sig./Nbkg_qU_sig;

fprintf('--------------------------------\n');
fprintf('Tot (analysis interval, only signal points) signal to background   %.2f \n',Nsig/Nbkg_sig);
fprintf('Mean signal to background   %.2f \n',mean(SB_qU_sig));
fprintf('Min signal to background   %.2f \n',min(SB_qU_sig));
fprintf('Max signal to background   %.2f \n',max(SB_qU_sig));

% signal-to-background equilibrium
qU = R.RunData.qU(R.exclDataStart:end);
Sig2Bkg1_lin = interp1(SB_qU(14:20),qU(14:20)-R.ModelObj.Q_i,1,'lin'); %cross check with linear interp (more robust, less accurate)
Sig2Bkg1 = interp1(SB_qU(14:18),qU(14:18)-R.ModelObj.Q_i,1,'spline');
fprintf('Signal to background  equilibrium  at qU-E0 = %.2f eV\n',Sig2Bkg1);


%% plot

f1 = figure('Units','normalized','Position',[1.1,1.1,0.5,0.4]);%  GetFigure;
bar(qU-R.ModelObj.Q,SB_qU,'FaceColor',rgb('Silver'),'EdgeColor',rgb('Silver'));
hold on;
p1 = plot(linspace(-40,0,27),ones(27,1),':k','LineWidth',2.5);
xlim([-40.5 0]);
ylim([0 700]);
xlabel(sprintf('Retarding energy - {\\itE}_0 (eV)'));
ylabel('Signal-to-background ratio');
set(gca,'YScale','log');

leg = legend(p1,sprintf('Signal-to-background equilibrium at {\\itqU - E}_0 = %.0f eV',Sig2Bkg1),'Location','northeast');
legend boxoff
PrettyFigureFormat('FontSize',22);
leg.FontSize = get(gca,'FontSize');

ax = gca;
ax.XAxis.Exponent = 0;
% %%
pltdir = [getenv('SamakPath'),'knm2ana/knm2_background/plots/'];
pltfile = sprintf('%sknm2_Sig2Bkg.pdf',pltdir);
export_fig(pltfile);


%% PRD review - reproduce
%Counts_qU_AnaRange = R.RunData.TBDIS(R.exclDataStart:end);
%Counts_AnaRange = sum(Counts_qU_AnaRange);

Idx = 1:numel(R.RunData.qU); %R.exclDataStart
t_i = (R.RunData.qUfrac(Idx).*R.RunData.TimeSec);
R_i = R.RunData.TBDIS(Idx)./t_i;%R.exclDataStart:end
R_i_bkg = 0.22.*ones(numel(R_i),1); %cps
R_bkg = sum(R_i_bkg);
R_i_sig = R_i-R_i_bkg;
R_sig = sum(R_i_sig);

%calculate signal to background with summed rates (what was done in PRD review)
Sig2Bkg_Rates = sum(R_i_sig)./mean(R_i_bkg);

%calculate signal to background with summed counts
Sig2Bkg_Counts = sum(R_i_sig.*t_i)./sum(R_i_bkg.*t_i);
Sig2Bkg_Counts_qU = (R_i_sig.*t_i)./(R_i_bkg.*t_i);

%% interpolate 
qU_all = R.RunData.qU-R.ModelObj.Q;
qU_inter = [-90:10:-10];
qU_all_inter = linspace(-90,-10,1e3);
Sig2Bkg_Counts_qU_all_inter = interp1(qU_all,Sig2Bkg_Counts_qU,qU_all_inter,'spline');
Sig2Bkg_Counts_qU_inter = interp1(qU_all,Sig2Bkg_Counts_qU,qU_inter,'spline');
%
fprintf('Signal to background \n');
IdxDisp = [1,5,7,9,11,16,21,27];
for i=1:numel(qU_inter)%IdxDisp)
 fprintf('qU-E0 = %.0f eV, S/B = %.1f (inter)\n',qU_inter(i),Sig2Bkg_Counts_qU_inter(i));
%fprintf('qU-18574 = %.0f eV, S/B = %.1f \n',qU_all(IdxDisp(i)),Sig2Bkg_Counts_qU(IdxDisp(i)));
end


%% full intervall (-300)

t_i = [R.RunData.qUfrac;R.RunData.qUfrac_RM].*R.RunData.TimeSec;
R_i = [R.RunData.TBDIS;R.RunData.TBDIS_RM]./t_i;%R.exclDataStart:end
R_i_bkg = 0.22.*ones(numel(R_i),1); %cps
R_bkg = sum(R_i_bkg);
R_i_sig = R_i-R_i_bkg;
R_sig = sum(R_i_sig);
Sig2Bkg_Counts = sum(R_i_sig.*t_i)./sum(R_i_bkg.*t_i);



%% 
SpectrumCounts = R.RunData.TBDIS(R.exclDataStart:end);
BkgCounts = 0.22.*(R.RunData.qUfrac(R.exclDataStart:end).*R.RunData.TimeSec); %cps
SignalCounts = SpectrumCounts-BkgCounts;
SignalCounts_sum = sum(SignalCounts);
Sig2Bkg_Counts = SignalCounts_sum/sum(BkgCounts);







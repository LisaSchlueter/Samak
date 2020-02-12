%
% Compute a Covariance Matrix 
% Apply to Krypton 83m Analysis 
% Moosfitter data challenge
%
% Th. Lasserre - CEA Saclay
% January 2018
%

format long;

nSpectra = 1e3;
norm = 10000;

global KrObjSim;
KrObjSim =krl332minuit_init('DopplerEffectFlag','OFF');
KrObjSim.HVRipples     = 'ON';
KrObjSim.HVRipplesP2PV = 0.52;
KrObjSim.TimeYear      = 1/365.25/86400*norm ; KrObjSim.TimeSec = 1*norm;
KrObjSim.qUfrac        = ones(numel(KrObjSim.qUfrac),1);
KrObjSim.L3_32_E_i     = 30472.334;
KrObjSim.L3_32_W_i     = 1.033167 ;
KrObjSim.L3_32_Phi0_i  = 726.342  ;
KrObjSim.BKG_RateSec_i = 276.413  ;
KrObjSim.ComputeKrDS; KrObjSim.ComputeKrIS;
KrObjSim.DisplayKrInfo();

% plot initial spectra
figure(1)
hold on
plot(KrObjSim.qU, KrObjSim.KrIS,'Color','Black','LineWidth',2)
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');

HVRipplesValue = ...
      KrObjSim.HVRipplesP2PV ...
    + KrObjSim.HVRipplesP2PV/20*randn(nSpectra,1); % Scan HVRipple Values
myKrIS      = zeros(nSpectra,KrObjSim.nqU);        % for storing IS and ISE
myKrISE     = zeros(nSpectra,KrObjSim.nqU);

% hv ripple distribution
figure(2)
nhist(HVRipplesValue)
xlabel('HV ripple amplitude (V)','FontSize',14);

% loop
tic
progressbar('Generate KrIS');
for i= 1:nSpectra       
progressbar(i/nSpectra);

% change hv ripple and compute integral spectrum
KrObjSim.HVRipplesP2PV = HVRipplesValue(i);
KrObjSim.ComputeKrIS;
myKrIS(i,:)  = KrObjSim.KrIS;
myKrISE(i,:) = KrObjSim.KrISE;
end
hold off

toc

% compute covariance matrix 
MyCovMat = cov(myKrIS);

%% plot covariance matrix
figure(3)
subplot(2,2,1)
imagesc(MyCovMat);
colormap(copper)
colorbar
title('HV Ripple CM - No Stat. Fluct.')

subplot(2,2,2)
surf(MyCovMat);
colormap(copper)
colorbar
title('HV Ripple CM - No Stat. Fluct.')

subplot(2,2,3)
SE2 = diag(mean(myKrIS(:,:)));
imagesc((SE2));
colormap(copper)
colorbar
title('Stat. Fluct. CM')

subplot(2,2,4)
M = (SE2) + MyCovMat;
imagesc(M);
colormap(copper)
colorbar
title('HV Ripple CM - With Stat. Fluct.')

%% Fractional Covariance Matrix 
SS = sqrt(diag(MyCovMat));
MyCovMatFrac = MyCovMat./(SS*SS');

%% plot spectra
figure(4)
scaleerr = 1;
StatE = sqrt(mean(myKrIS(:,:)))*scaleerr;
SystE = sqrt(diag(MyCovMat)')*scaleerr;
TotE  = sqrt(StatE.^2+SystE.^2);
hold on
hstat = errorbar(KrObjSim.qU,mean(myKrIS(:,:)),StatE,'Color','Black')
hold on
hsyst = errorbar(KrObjSim.qU,mean(myKrIS(:,:)),SystE,'Color','Red');
htot  = errorbar(KrObjSim.qU,mean(myKrIS(:,:)),TotE,'Color','Blue');
hold off
legend([hstat hsyst htot],'Stat. Error','HVR Syst. Error','Tot. Error','Location','NorthEast') ;
grid on

% figure(999)
% subplot(1,3,1)
% imagesc(MyCovMat);
% colorbar;
% subplot(1,3,2)
% imagesc(MyCovMatFrac.*(SS*SS'));
% colorbar;
% subplot(1,3,3)
% SS = sqrt(diag(MyCovMat))/100;
% imagesc(MyCovMatFrac.*(SS*SS'));
% colormap(copper);
% colorbar;

% Save Covariance Matrix
mystr = sprintf('./krl332minuit_fitcovm.mat');
save(mystr,'MyCovMat','MyCovMatFrac','-mat');
            
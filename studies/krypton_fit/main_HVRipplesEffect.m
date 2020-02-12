addpath(genpath('../../../Samak2.0'));
clear all;

%% Fit 
myPixel = 1;
generalFitOpt = {'CPS','ON', 'FPD_Segmentation','PIXEL','Pixel',myPixel, 'Fitter', 'Minuit','TD','KrL3_32_HS',...
           'DopplerEffectFlag','ON','ConvFlag','DopplerEffect', 'HVRipples','ON','mypub','OFF'};
nfit = 100;  %nmb of fits
npar = 4;     %nmb of parameters
mypar = zeros(nfit,npar+1);   %Matrix to store parameters
myerr = zeros(nfit,npar+1);   %Uncertainties
mychi2min = zeros(nfit, 1); %Chi2
HVRipplesValue = [0.01:0.01:1];%Scan HVRipple Values in 0.01 V steps

tic
progressbar('Krypton Line Fits');
for i = 1:nfit
progressbar(i/nfit);
[par err chi2min ndof] = KrLineFitData(generalFitOpt{:},'HVRipplesP2PV', HVRipplesValue(i));  %perform fit with different HVRipplesValue (Model)
mypar(i,:)=par(:);
myerr(i,:)=err(:);
mychi2min(i,:)=chi2min;
end
toc

%% Plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(100);
plot(HVRipplesValue*1e03,mychi2min,'-');%./ndof,'-');
xlabel('HV-Ripple (mV)','FontSize',16);
ylabel('$\mathbf{\chi^2}$','Interpreter','latex','FontSize',16); % /\textrm{\textbf{ndof}}
myleg = sprintf('Fit to 83m-Kr L3-Line \n Pixel %i', myPixel);
l1 = legend(myleg, 'Location', 'northwest');
l1.FontSize =16;
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A; %Object of Kr class, used for fit
 FitE       =   A.L3_32_E_i+mypar(:,1);
 FitEerr    =   myerr(:,1);
 FitW       =   A.L3_32_W_i+mypar(:,2)*1e3;
 FitWerr    =   myerr(:,2)*1e3;
 FitPhi0    =   A.L3_32_Phi0_i+mypar(:,3);
 FitPhi0err =   myerr(:,3);
 FitOff     =  (A.BKG_RateSec_i+mypar(:,4))*(A.TimeSec.*A.qUfrac(1));
 FitOfferr  = myerr(:,4)/(A.TimeSec.*A.qUfrac(1));

% figure(200);
% %title('Fit to 83m-Kr L3-Line Pixel 1');
% subplot(2,2,1);
% hist(FitE*1e-03);
% xlabel('Fit: Line Energy');
% subplot(2,2,2);
% hist(FitW);
% xlabel('Fit: Line Width');
% subplot(2,2,3);
% hist(FitPhi0);
% xlabel('Fit: Line Activity');
% subplot(2,2,4);
% hist(FitOff);
% xlabel('Fit: Offset');
% 
% figure(300);
% %title('Fit to 83m-Kr L3-Line Pixel 1');
% subplot(2,2,1);
% plot(HVRipplesValue*1e03, FitE*1e-03);
% xlabel('HV-Ripple (mV)','FontSize',16);
% ylabel('Energy (keV)');
% grid on;
% subplot(2,2,2);
% plot(HVRipplesValue*1e03, FitW)
% xlabel('HV-Ripple (mV)','FontSize',16);
% ylabel('Width (meV)');
% grid on;
% subplot(2,2,3);
% plot(HVRipplesValue*1e03, FitPhi0)
% xlabel('HV-Ripple (mV)','FontSize',16);
% ylabel('Activity (cps)');
% grid on;
% subplot(2,2,4);
% plot(HVRipplesValue*1e03, FitOff);
% xlabel('HV-Ripple (mV)','FontSize',16);
% ylabel('Offset (cps)');
% grid on;
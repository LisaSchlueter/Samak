addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];
myFitOpt = {'Mode', 'Data', 'TD', 'KrL3_32_Satellites',...
 'Fitter','Minuit', 'FPD_Segmentation','PIXEL', ...
'Pixel',1 ,'DopplerEffectFlag','OFF', 'HVRipples','OFF',...
'ConvFlag','DopplerEffect', 'MultiPeaksFlag','ON'};

%Fit results
[mypar myerr mychi2min myndof] = KrLineFit_L3_32_Satellites(myFitOpt{:});
FitKrIS = KrLineModel7par_sat(mypar); %KrIS with parameter from fit

%Data
%[qU krdatamap] = readKrData('TD',TD);
filename= '../../krypton-data/33196_L3-32_high_stat_satellite/33196_channel0.txt';
mydata_tmp = importdata(filename);
mydata(:,1) = flip(mydata_tmp(:,1));
mydata(:,2) = flip(mydata_tmp(:,2));
mydata(:,3) = flip(mydata_tmp(:,3));

figure(1);%satelitte FIT

%% 
set(gcf,'Color','w') %background white
subpl1 = subplot(2,1,1);
errorbar(mydata(:,1)*1e-3, mydata(:,2), mydata(:,3), 'x red');
hold on;
plot(mydata(:,1)*1e-3, FitKrIS,'-','Color', orange);
xlabel('qU [keV]');
ylabel('N');
xlim([min(mydata(:,1)*1e-3) max(mydata(:,1)*1e-3)]);
grid on;
hold off;
legend('Data L3-32-HS (Pixel 1)','Fit','Location', 'southwest');
title('Krypton L3-32-HS Pixel 1');

subpl2 = subplot(2,1,2);
errorbar(mydata(:,1)*1e-3, mydata(:,2)-FitKrIS, mydata(:,3),'x', 'Color', orange);
hold on;
myRefLine= zeros(size(mydata(:,1),1), 1);
plot(mydata(:,1)*1e-3, myRefLine, '--k');
xlabel('qU [keV]');
ylabel('Residuals');
xlim([min(mydata(:,1)*1e-3) max(mydata(:,1)*1e-3)]);
grid on; 
hold off;
legend('Data-Fit');
%export_fig ../krypton_fit_LisaTest/plots/Fit_L3-32_HS_new.pdf



% subplot(3,1,3); %figure with 3 plots
% errorbar(mydata(:,1)*1e-3, (mydata(:,2)-FitKrIS)./mydata(:,2), (FitKrIS./(mydata(:,2).^2)).*mydata(:,3),'x', 'Color', orange);
% hold on;
% myRefLine= zeros(size(mydata(:,1),1), 1);
% plot(mydata(:,1)*1e-3, myRefLine, '--k');
% xlabel('qU [keV]');ylabel('Residuals/Data');
% xlim([min(mydata(:,1)*1e-3) max(mydata(:,1)*1e-3)]);
% grid on; 
% hold off;
% legend('L3-32-HS Pixel 1');
 

% % figure(2); %relative residuals
% errorbar(mydata(:,1)*1e-3, (mydata(:,2)-FitKrIS)./mydata(:,2), (FitKrIS./(mydata(:,2).^2)).*mydata(:,3),'x', 'Color', orange);
% hold on;
% myRefLine= zeros(size(mydata(:,1),1), 1);
% plot(mydata(:,1)*1e-3, myRefLine, '--k');
% xlabel('qU [keV]');ylabel('Residuals/Data');
% xlim([min(mydata(:,1)*1e-3) max(mydata(:,1)*1e-3)]);
% grid on; 
% hold off;
% legend('L3-32-HS Pixel 1');
% export_fig /home/lisa/Dokumente/Studium/Masterstudium/Masterarbeit/KATRIN_MatlabCode/plots/Krypton_firstSim/Fit_L3-32_HS_residuals_new.pdf

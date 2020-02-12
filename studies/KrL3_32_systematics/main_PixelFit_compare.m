addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];
clear all;
format long;

%Load Fit Results
load('FitResults_Pixels_NoCovMat.mat');
load('FitResults_Pixels_CovMat.mat');

ndof = A{1,1}.nqU-4; 
        
% Exclude Pixel
Pixel = [1:148];  
    parfor i= 1:148 
 if max(krdatamap(i,1,:))<10  || chi2pvalue(mychi2min(i), ndof) < 0.05 || chi2pvalue(mychi2min_NoCovMat(i), ndof) < 0.05
    fprintf('Pixel excluded: %u \n', Pixel(i))
    E_Fit(i) = NaN;
    E_FitError(i) = NaN;
    W_Fit(i) = NaN;
    W_FitError(i) = NaN;
    Phi0_Fit(i) = NaN;
    Phi0_FitError(i) = NaN;
    Offset_Fit(i) = NaN;
    Offset_FitError(i) = NaN;
    Pixel(i) = NaN; 
    mychi2min(i) = NaN; 
    
    E_Fit_NoCovMat(i) = NaN;
    E_FitError_NoCovMat(i) = NaN;
    W_Fit_NoCovMat(i) = NaN;
    W_FitError_NoCovMat(i) = NaN;
    Phi0_Fit_NoCovMat(i) = NaN;
    Phi0_FitError_NoCovMat(i) = NaN;
    Offset_Fit_NoCovMat(i) = NaN;
    Offset_FitError_NoCovMat(i) = NaN;
    mychi2min_NoCovMat(i) = NaN;
 end
    end 

%% Plots 

% Effect on Width and Energy 
f = figure(1001); 
set(f, 'Units', 'normalized', 'Position', [10*0.2, 10*0.1, 10*0.9, 10*0.8]);  
%set(gcf,'Color','w') %background white
subplot(2,2,1);
h10 = histogram(W_Fit-W_Fit_NoCovMat);
mygreen = rgb('DarkGreen'); h10.FaceColor = mygreen;
xlabel('$\Delta$ \textrm{\textbf{Width (meV)}}', 'Interpreter', 'latex', 'FontSize', 18);
legend(sprintf('CovMat \n - NoCovMat'));
%grid on;

subplot(2,2,2);
h1 = histogram(W_FitError-W_FitError_NoCovMat);
mygreen = rgb('DarkGreen'); h1.FaceColor = mygreen;
xlabel('$\Delta$ \textrm{\textbf{Width Uncertainty (meV)}}', 'Interpreter', 'latex','FontSize', 18);
%grid on;
%legend('CovMat - NoCovMat');

subplot(2,2,3);
h1 = histogram(E_Fit-E_Fit_NoCovMat);
mygreen = rgb('DarkGreen'); h1.FaceColor = mygreen;
xlabel('$\Delta$ \textrm{\textbf{Energy (eV)}}', 'Interpreter', 'latex','FontSize', 18);
%grid on;
%legend('CovMat - NoCovMat');

subplot(2,2,4);
h1 = histogram(E_FitError-E_FitError_NoCovMat);
mygreen = rgb('DarkGreen'); h1.FaceColor = mygreen;
xlabel('$\Delta$ \textrm{\textbf{Energy Uncertainty (eV)}}', 'Interpreter', 'latex','FontSize', 18);
%legend('CovMat - NoCovMat');
%grid on;
%export_fig './plots/Comparison/fig/Comparison_WidthEnergy.fig'
export_fig './plots/Comparison/Comparison_WidthEnergy.png'

%%%%%  Systematic Error on Width
f2 = figure(234); 
set(f2, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.8, 0.8]);  
set(gcf,'Color','w') %background white
WsysError = sqrt(W_FitError.^2-W_FitError_NoCovMat.^2);
plot(Pixel, WsysError,'x');
xlabel('Pixel');
ylabel('$\mathbf{\sigma}_{\textrm{\textbf{sys}}} \textrm{\textbf{(meV)}}$', 'Interpreter', 'latex', 'FontSize', 16);
xticks([1:9:145]);
title(sprintf('KATRIN Gaseous Krypton 83m L3-32 Line \n Systematic Uncertainty on Width due to Gaussian HV-Ripple Variation of 50%'));
grid on;
 %export_fig './plots/Comparison/fig/WidthSysUncertaintyRipple.fig'
 %export_fig './plots/Comparison/WidthSysUncertaintyRipple.pdf'

%%% Systematic Error on Phi0 and Offset
OffsetsysError = sqrt(Offset_FitError.^2-Offset_FitError_NoCovMat.^2);
Phi0sysError = sqrt(Phi0_FitError.^2-Phi0_FitError_NoCovMat.^2);  
plot(Pixel, Phi0sysError,'x');
xlabel('Pixel');
ylabel('$\mathbf{\sigma}_{\textrm{\textbf{sys}}} \textrm{\textbf{(cps)}}$', 'Interpreter', 'latex', 'FontSize', 16);
xticks([1:9:145]);
title(sprintf('KATRIN Gaseous Krypton 83m L3-32 Line \n Systematic Uncertainty on Activity due to Gaussian HV-Ripple Variation of 50%'));
grid on;
% export_fig './plots/Comparison/Phi0SysUncertaintyRipple.pdf'
% export_fig './plots/Comparison/fig/Phi0SysUncertaintyRipple.fig'

%FPDViewer((W_FitError-W_FitError_NoCovMat));
%FPDViewer((W_Fit-W_Fit_NoCovMat));


% Effect on Phi0 and Offset
f = figure(1002); 
set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.8, 0.8]);  
set(gcf,'Color','w') %background white
subplot(2,2,1);
h10 = histogram(Offset_Fit-Offset_Fit_NoCovMat);
xlabel('$\Delta$ \textrm{\textbf{Offset (cps)}}', 'Interpreter', 'latex');
legend('CovMat - NoCovMat');
%grid on;

subplot(2,2,2);
h1 = histogram(Offset_FitError-Offset_FitError_NoCovMat);
xlabel('$\Delta$ \textrm{\textbf{Offset Uncertainty (cps)}}', 'Interpreter', 'latex');
%grid on;
legend('CovMat - NoCovMat');

subplot(2,2,3);
h1 = histogram(Phi0_Fit-Phi0_Fit_NoCovMat);
xlabel('$\Delta$ \textrm{\textbf{Phi0 (cps)}}', 'Interpreter', 'latex');
%grid on;
legend('CovMat - NoCovMat');

subplot(2,2,4);
h1 = histogram(Phi0_FitError-Phi0_FitError_NoCovMat);
xlabel('$\Delta$ \textrm{\textbf{Phi0 Uncertainty (cps)}}', 'Interpreter', 'latex');
legend('CovMat - NoCovMat');
%grid on;
% export_fig './plots/Comparison/fig/Comparison_OffsetPhi0.fig'
% export_fig './plots/Comparison/Comparison_OffsetPhi0.png'
                
%%%%%%Plot Fit Result: absolute values

% figure(1)
% subplot(2,2,1);
% h10 = histogram(W_Fit_NoCovMat)
% h1.FaceColor = 'b';
% hold on;
% h2= histogram(W_Fit); 
% h2.FaceColor = 'c';
% xlabel('Width (meV)');
% grid on;
% legend('NoCovMat', 'CovMat');
% hold off;
% 
% subplot(2,2,2);
% h1 = histogram(W_FitError_NoCovMat)
% h1.FaceColor = 'b';
% hold on;
% h2= histogram(W_FitError); 
% h2.FaceColor = 'c';
% xlabel('Width Uncertainty(meV)');
% grid on;
% legend('NoCovMat', 'CovMat');
% hold off;
% 
% subplot(2,2,3);
% h1 = histogram(E_Fit-E_Fit_NoCovMat)
% h1.FaceColor = 'b';
% hold on;
% h2= histogram(E_Fit); 
% h2.FaceColor = 'c';
% xlabel('Energy (eV)');
% grid on;
% legend('NoCovMat','CovMat');
% hold off;
% 
% subplot(2,2,4);
% h1 = histogram(E_FitError_NoCovMat)
% h1.FaceColor = 'b';
% hold on;
% h2= histogram(E_FitError); 
% h2.FaceColor = 'c';
% xlabel('Energy Uncertainty(eV)');
% grid on;
% legend('NoCovMat', 'CovMat');
% hold off;


%%%%%%% FPDViewer
%figure(99);
%errorbar(Pixel,E_Fit, E_FitError,'x'); hold on;
%errorbar(Pixel, E_Fit_NoCovMat, E_FitError_NoCovMat,'s'); hold off;
% grid on;
% xlabel('Pixel');
% ylabel('Energy (eV)');

% figure(11);
% h1 = histogram(myfit_xi2min_NoCovMat./ndof);
% hold on;
% h2 = histogram(mychi2min./ndof); hold off;
    

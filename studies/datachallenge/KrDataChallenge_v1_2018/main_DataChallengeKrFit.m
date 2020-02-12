%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of Samak Spectrum for KrDataChallenge
% Fit of other Spectra
% Comparison 
%
% Lisa Schlueter 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];
Modelname= 'Samak';
% Creation of a integral spectrum with my "true parameters" + fluctuations
% (has to be done only once) 
createSpectra = 'OFF';
switch createSpectra 
    case 'ON'
Obj_KrDC=InitKATRIN_KrDataChallenge_benchmark();
Obj_KrDC.FPD_Pixel = 1;
Obj_KrDC.L3_32_W_i = 1.19; %Cannot init this in InitKrKATRIN
Obj_KrDC.L3_32_E_i = 30472.5;
Obj_KrDC.ComputeKrDS;
Obj_KrDC.ComputeKrIS;
Obj_KrDC.AddStatFluctKrIS(); %my true parameters + stat fluctuations
SamakSpec = [Obj_KrDC.qU , Obj_KrDC.KrIS , Obj_KrDC.KrISE];
format long;
save('../../krypton-data/DataChallenge18V1/kdc_v1_samak.txt','SamakSpec','-ascii');
end 

%%%%%%%%% Fit Spectra (created by Samak, SCC, Martin,....) %%%%%%%%%%%%%%
A=InitKATRIN_KrDataChallenge_benchmark(); %benchmark parameter setting
A.L3_32_W_i = 1.19; %Cannot init this in InitKrKATRIN
A.L3_32_E_i = 30472.5;
[par err chi2min] = KrL3LineFitData_NoCovMat('KrObject', A,'plots', 'OFF'); %still need to change ModelName in readKrData!
ndof = A.nqU-4;

% Reading Data for plot
[qU, krdatamap] = readKrData('TD','DummyKrdc1'); 
Count = krdatamap(1,1,:);
CountErr = krdatamap(1,2,:);

% Plot Fit Results
E_Fit             = (A.L3_32_E_i+par(1));      %[eV]
E_FitError        = (err(1));
W_Fit             = (A.L3_32_W_i+par(2))*1e3;  %[meV]
W_FitError        = (err(2))*1e3;
Phi0_Fit          = (A.L3_32_Phi0_i+par(3));   %[cps]
Phi0_FitError     = err(3);
Offset_Fit        = (A.BKG_RateSec_i+par(4));  %[cps]
Offset_FitError   = err(4);

%True Parameter
% A.L3_32_W_i = 1.22;
% A.L3_32_E_i = 30472.7;
% A.L3_32_Phi0_i=55;
% A.BKG_RateSec_i=7;
%par = zeros(5,1);

fitleg = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f\n E= %.3f ± %.3f eV \n W=%.0f  ± %.0f meV \n A=%.3f  ± %.3f cps \n O=%.3f  ± %.3f cps',...
    chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError, Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
fig = figure('NumberTitle','off','rend','painters' ,'pos',[10 10 1300 1000]);
%set(gcf,'Color','w') %background white

s1= subplot(2,1,1); %Data + Fit
%gdata = errorbar(qU(:)*1e-03, Count(:),CountErr(:), 's','MarkerSize',2, 'Color', orange);
ModelErr = sqrt(KrLineModel4par(par,A).*(A.TimeSec.*A.qUfrac(1)))./(A.TimeSec.*A.qUfrac(1));
gdata = errorbar(qU(:)*1e-03, Count(:),ModelErr, 's','MarkerSize',2, 'Color', orange);
hold on;
gfit = plot(qU(:)*1e-03, KrLineModel4par(par,A),'-', 'Color', 'Red');
%gfit = errorbar(qU(:)*1e-03, KrLineModel4par(par,A),ModelErr,'x', 'Color', 'Red');
offline = line([qU(1)*1e-03, qU(end)*1e-03], [Offset_Fit, Offset_Fit],'LineStyle','--', 'LineWidth', 1.5);
hold off;
xlim([30.4670 30.475]);
xlabel('qU (keV)','FontSize', 14);
ylabel('$\dot{\textrm{\textbf{N}}}$ \textbf{(cps)}','Interpreter','latex','FontSize', 14);
dataleg = sprintf('Data (%s)', Modelname);
intleg = legend([gdata gfit offline], dataleg, fitleg, 'Offset');
%title(sprintf('Krypton Data Challenge 2018 V1: \n SSC (No stat. fluctuations; 1 day long subruns) and Samak Model (All Parameters fixed)'));
grid on;

s2= subplot(2,1,2); %Residuals
plot(qU(:)*1e-03, (Count(:)-KrLineModel4par(par,A))./KrLineModel4par(par,A),'x' ,'Color', orange );
%errorbar(qU(:)*1e-03, (Count(:)-KrLineModel4par(par), CountErr(:), 'x', 'Color', orange )
hold on;
refline = line([qU(1)*1e-03, qU(end)*1e-03], [0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
hold off;
grid on;
xlim([30.4670 30.475]);
ylabel('\textbf{Residuals}/Counts', 'FontSize',14, 'Interpreter', 'latex');
xlabel('qU (keV)','FontSize', 14);
export_fig ./plots/fig/KrDataChallenge_samak.fig
export_fig ./plots/KrDataChallenge_samak.pdf

%%  Short study of new relativistic TF
%  figure(001);
%  set(gcf,'Color','w') %background white
% % %subp1 = subplot(2,1,1);
%  qUtest = 3.04676e4;
% % %global A;
% myTe = [30467.0:0.1:30475.4];
%  plot(myTe-qUtest, A.KTF(myTe,qUtest,A.MACE_R_eV),'x');
% % hold on;
% % plot(Obj_KrDC2.Te-qUtest, Obj_KrDC2.KTF(Obj_KrDC2.Te,qUtest,Obj_KrDC2.MACE_R_eV),'--');
% % hold off;
% % leg1 = sprintf('Samak TF: \n relativitic');
% % legend(leg1, 'non relativistic', 'Location', 'NorthWest');
% % title('Relativistic Transmission Function (Samak)');
% % legend(sprintf('qU_0 = %u eV',qUtest), 'Location', 'southeast');
%  grid on;
% % xlim([-0.5 3]);
%  xlabel('$$\textrm{\textbf{T}}_{\textrm{\textbf{e}}}-\textrm{\textbf{q}}\textrm{\textbf{U}}_0$$ (\textrm{\textbf{eV}})' , 'Interpreter', 'latex','FontSize', 16);
%  ylabel('\textbf{KTF(T}$$_{\textrm{\textbf{e}}}$$, \textbf{qU}$$_0$$)', 'Interpreter', 'latex');
% % ylim([-0.1 1.1]);
% 
%  SamakKTF = [(myTe-qUtest) ;  A.KTF(myTe,qUtest,A.MACE_R_eV)]';
%  save('Samak_relKTF.txt', 'SamakKTF', '-ascii');
%  export_fig 'KrDataChallenge_v1_2018/Samak_relKTF.png'

addpath(genpath('../../../Samak2.0'));

myFitOpt = {'Mode', 'Sim', 'TD', 'MoSFitterJuly2017', 'CovMat', 'OFF', 'MultiPeaksFlag', 'OFF','FPD_Segmentation','OFF','HVRipples','OFF','ConvFlag','OFF','CPS','OFF'};

[mypar myerr mychi2min myndof] = KrLineFit_Lisa(myFitOpt{:});
global A; %defined in KrLineFit.m: A ist object of class Kr with "true"+stat. fluctuation values (in case of Sim; otherwise true values)

myObj =InitKrKATRIN_LisaFit('TD','MoSFitterJuly2017','FPD_Segmentation','OFF','HVRipples','OFF','ConvFlag','OFF','CPS','OFF','MultiPeaksFlag', 'OFF'); %"true values"

figure(1);
myObj.ComputeKrDS(); myObj.ComputeKrIS();
ModelKrIS = myObj.KrIS(); %true values no fluctuation 
ModelKrDS = myObj.KrDS();
%ModelKrIS = A.KrIS; %true values + stat. Fluctation

myObj.ComputeKrDS(...
    'E_bias',mypar(1),...
    'W_bias',mypar(2),...
    'Phi0_bias',mypar(3),...
    'B_bias',mypar(4)); 
myObj.ComputeKrIS();
myFitKrIS= myObj.KrIS(); %fit results

subplot(2,1,1);
title('Comparison Simulation - Fit');
plot(myObj.qU*1e-03, ModelKrIS, 'x');
hold on;
plot(myObj.qU*1e-03, myFitKrIS);
hold off;
grid on;
ylabel('$\dot{\textrm{\textbf{N}}}$', 'Interpreter', 'latex');
xlim([min(myObj.qU)*1e-03 max(myObj.qU)*1e-03]);
xlabel('qU [keV]');
legend('Simulation','Fit to (Simulation + stat. Fluctuations)')
%myFitKrIS= KrLineModel4par(mypar);
subplot(2,1,2);
plot (myObj.qU*1e-03, (ModelKrIS-myFitKrIS),'x');
grid on;
%ylim([-0.000001 0.000001]);
xlim([min(myObj.qU)*1e-03 max(myObj.qU)*1e-03]);
ylabel('Residuals');
xlabel('qU [keV]');
%export_fig ../krypton_fit_LisaTest/plots/Fit_Simulation_new.pdf
%
figure(2);
plot(myObj.Te*1e-03, ModelKrDS);





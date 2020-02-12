addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];
myFitOpt = {'Mode', 'MoSFitterJuly2017', 'TD', 'MoSFitterJuly2017',...
 'CovMat', 'OFF', 'Fitter','Minuit', 'FPD_Segmentation','OFF', ...
'Pixel',1 ,'DopplerEffectFlag','OFF', 'HVRipples','OFF', 'MultiPeaksFlag','OFF'};

%Fit results
[mypar myerr mychi2min myndof] = KrLineFit_Lisa(myFitOpt{:});
FitKrIS = KrLineModel4par(mypar); %KrIS with parameter from fit

%Data MosFItter
filename = '../../krypton-data/MoSFitterJuly2017/output_int_L3-32.txt';

mydata = importdata(filename);

%Plot
subpl1 = subplot(2,1,1);
errorbar(mydata(:,1)*1e-3, mydata(:,2), mydata(:,3), 'x red');
hold on;
plot(mydata(:,1)*1e-3, FitKrIS,'-','Color', orange);
xlabel('qU [keV]');
ylabel('$\dot{\textrm{\textbf{N}}}$', 'Interpreter', 'latex');
xlim([min(mydata(:,1)*1e-3) max(mydata(:,1)*1e-3)]);
grid on;
hold off;
dim = [.65 .6 .3 .3];
str = {'\textbf{Fit Result:}','A = 84.5 cps','E0 = 30472.3 eV','W = 1314 meV', 'O = 32.0 cps'};
t= annotation('textbox', dim, 'String', str,'FitBoxToText','on');
t.FontSize =14;
t.HorizontalAlignment = 'center';
t.Interpreter = 'latex'
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
%linkaxes([subpl1,subpl2],'x');
%export_fig ../krypton_fit_LisaTest/plots/Fit_MoSFitterJuly2017_new.pdf


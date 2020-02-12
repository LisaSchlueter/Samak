% 
% Divide the FPD detector within four rings
% Perform the final fit on each ring
% Display fit Results for Each Ring
% Analysis Perfomed
% - With Stat Only
% - With Stat+Sys No correlations between rings
% 
% T. Lasserre
% Last Updated: 12/07/2019
%
RunList               = 'KNM1';
exclDataStart         = 2;
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
BkgCM                 = 'ON';

%% First Step, Determine the common endpoint for the Uniform Fit
MC = MultiRunAnalysis('RunList','KNM1',...
    'RingMerge','Full',...
    'fixPar','5 6 7 8 9 10 11',...
    'exclDataStart',exclDataStart,...
    'DataType','Real',...
    'chi2','chi2CMShape');
MC.ComputeCM('SysEffects',SysEffects,'BkgCM','ON');
MC.Fit;

%% Second Step, fix Common endpoint and fit neutrino mass squared per ring
MC4 = MultiRunAnalysis('RunList','KNM1',...
    'RingMerge','Full',...
    'fixPar','2 5 6 7 8 9 10',...
    'exclDataStart',exclDataStart,...
    'DataType','Real',...
    'chi2','chi2CMShape',...
    'i_Q', MC.FitResult.par(2));
MC4.ComputeCM('SysEffects',SysEffects,'BkgCM','ON');
MC4.GetPlotColor
MC4.SimulateStackRuns;
RC4 = RingAnalysis('RunAnaObj',MC4,'RingList',1:4);
RC4.FitRings('SaveResult','OFF','RecomputeFlag','ON');

%% Plot All Results

% Plot m2
RC4.PlotFits('SavePlot','ON','PlotPar',1,'Blind','OFF','PlotMode','Abs');

% Plot E0
RC4.PlotFits('SavePlot','ON','PlotPar',2,'Blind','OFF','PlotMode','Abs');

% Plot B
RC4.PlotFits('SavePlot','ON','PlotPar',3,'Blind','OFF','PlotMode','Abs');

% Plot N
RC4.PlotFits('SavePlot','ON','PlotPar',4,'Blind','OFF','PlotMode','Abs');

% Plot Offsets
RC4.PlotFits('SavePlot','ON','PlotPar',11,'Blind','OFF','PlotMode','Abs');

%% Ellispses
fig111 = figure(111);
set(fig111, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.7]);

%Ring 1,2,3,4
for i=1:4
cov_m2e0 = RC4.MultiObj(i).FitResult.errmat(1:2,1:2);
fit_m2e0 = [RC4.MultiObj(i).FitResult.par(1) RC4.MultiObj(i).FitResult.par(2)];
hold on
CLsigma(2)='2';
r_ellipse = RC4.MultiObj(i).DrawEllipseError(fit_m2e0,cov_m2e0,'CL',CLsigma(2));
d(i) = plot(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),'-','LineWidth',4);
fill(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),rgb('IndianRed'),'FaceAlpha',0.1);
plot(fit_m2e0(1),fit_m2e0(2),'s','Color','Black','MarkerSize',8,'LineWidth',3);
end

% Uniform
cov_m2e0 = MC.FitResult.errmat(1:2,1:2);
fit_m2e0 = [MC.FitResult.par(1) MC.FitResult.par(2)];
hold on
CLsigma(2)='2';
r_ellipse = MC.DrawEllipseError(fit_m2e0,cov_m2e0,'CL',CLsigma(2));
du = plot(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),'--','LineWidth',4);
fill(r_ellipse(:,1) + fit_m2e0(1),r_ellipse(:,2) + fit_m2e0(2),rgb('IndianRed'),'FaceAlpha',0.1);
plot(fit_m2e0(1),fit_m2e0(2),'s','Color','Black','MarkerSize',8,'LineWidth',3);


hold off
strleg='Data';
strleg2='Fit (Gaussian Approximation) ';
leg = legend([du,d(1),d(2),d(3),d(4)],'Uniform','Pseudo-Ring 1 (2 \sigma)','Pseudo-Ring 2 (2 \sigma)','Pseudo-Ring 3 (2 \sigma)','Pseudo-Ring 4 (2 \sigma)','Location','NorthEast');
legend('boxoff');
xlabel('m^2 (eV)','FontSize',30);
ylabel(sprintf('E_0 - %.2f (eV)',RC4.MultiObj(i).ModelObj.Q_i),'FontSize',30);
grid on;
PrettyFigureFormat
set(gca,'FontSize',30);
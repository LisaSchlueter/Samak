% 
% First Tritium Background Slope Constraint
% Th. Lasserre, April 2020
%

% start fit
starfit = 1;

% Import Data
FT_bg_knm2golden = importdata('./dataAnna/FT_bg_knm2golden.txt');

nbin = numel(FT_bg_knm2golden(starfit:end,1));
hv    = FT_bg_knm2golden(starfit:end,1);
c     = FT_bg_knm2golden(starfit:end,2);
t     = FT_bg_knm2golden(starfit:end,3);
% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
Data = [hv/1e3 c./t*1e3 sqrt(c)./t*1e3];

% Linear Model
parnames = {'baseline','slope'}; 
ParIni   = [mean(Data(:,2)) 0]; 

% Case 1) Fit Overall Range
% Fit
Args = {ParIni, Data, '-c', 'min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope    = %.3f ± %.3f\n',par(2),err(2));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));
fprintf('================================================================\n');

% Figure 1
close all
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 500])
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(:,1),Model(par,Data(:,1)),'LineWidth',3);
hfiteu =  plot(Data(:,1),Model(par-[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
hfited =  plot(Data(:,1),Model(par+[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
hold off;
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(1),err(1),par(2),err(2)),'Location','SouthWest')
grid on
publish_figurePDF(1,'FT_bg_knm2golden_1.pdf');

%% Case 2) Fit All Ranges
% Plot Slope & Error on Slope Verus HV-start
for i=1:1:(nbin-10)
DataHVscan = Data(i:end,:);
Args = {ParIni, DataHVscan, '-c', 'min; minos'};
[ parV(i,:), errV(i,:), chi2minV(i), errmatV(i,:,:) ] = fminuit('Chi2Gauss',Args{:});
HVmin(i) = DataHVscan(1,1);
end
% Figure 2
figure(2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
subplot(3,1,1)
hdata = errorbar(HVmin,parV(:,1),errV(:,1),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('baseline (mcps)','FontSize',14) 
grid on
PrettyFigureFormat;
subplot(3,1,2)
hdata = errorbar(HVmin,parV(:,2),errV(:,2),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('slope (mcps/keV)','FontSize',14) 
grid on
PrettyFigureFormat;
subplot(3,1,3)
hdata = semilogy(HVmin,parV(:,2)+1*errV(:,2),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('m+1\sigma (mcps/keV)','FontSize',14) 
grid on
PrettyFigureFormat;
publish_figurePDF(2,'FT_bg_knm2golden_2.pdf');

% Display Recommended Limit m+1sigma
SlopeLimit1sigma = table(HVmin',1*errV(:,2),norminv(0.68269,parV(:,2),errV(:,2)),'VariableNames',{'MinqU (keV)','CurrentLimit (mcps/keV)','RecommendedLimit (mcps/keV)'});
SlopeLimit1sigma = varfun(@(var) round(var, 2), SlopeLimit1sigma)

%% Fit Diagnostics
% Put the variables into a table, naming them appropriately
tbl = table(Data(:,1)-18.5737,Data(:,2),Data(:,3),'VariableNames',{'qUKEV','RateMCPS','RateErrorMCPS'});
% Specify and carry out the fit
mdl = fitlm(tbl,'RateMCPS ~ 1 + qUKEV','Weights',1./((Data(:,3)).^2),'RobustOpts','off','Exclude',[]);

% Leverage
figure(300)
plotDiagnostics(mdl,'leverage','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

% Cook Distance
figure(301)
plotDiagnostics(mdl,'cookd','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

% Dfbetas
figure(302)
plotDiagnostics(mdl,'Dfbetas','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

%% Change of Slope
figure(555)
set(gcf, 'Position',  [100, 100, 1000, 700])
subplot(2,1,1)
hdata = errorbar(Data(:,1)-18.5737,Data(:,2),Data(:,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('Rate (mcps)','FontSize',14) 
grid on
PrettyFigureFormat
subplot(2,1,2)
h=plot(Data(:,1)-18.5737,mdl.Diagnostics.Dfbetas(:,2) .* mdl.Coefficients.SE(2) ./ mdl.Coefficients.Estimate(2)*100,'ks',...
    'MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('\Delta Slope (%)','FontSize',14) 
grid on
legend([h],sprintf('all data slope \n %0.3f \\pm %0.3f mcps/keV',par(2),err(2)),'Location','SouthEast')
PrettyFigureFormat

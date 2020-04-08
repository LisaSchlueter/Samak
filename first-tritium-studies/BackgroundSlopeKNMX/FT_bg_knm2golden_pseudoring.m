% 
% First Tritium Background Slope Constraint
% Th. Lasserre, April 2020
%

% start fit
starfit = 1;

% Import Data
FT_bg_knm2golden_pseudoring0 = importdata('./dataAnna/FT_bg_knm2golden_pseudoring_0.txt');
FT_bg_knm2golden_pseudoring1 = importdata('./dataAnna/FT_bg_knm2golden_pseudoring_1.txt');
FT_bg_knm2golden_pseudoring2 = importdata('./dataAnna/FT_bg_knm2golden_pseudoring_2.txt');
FT_bg_knm2golden_pseudoring3 = importdata('./dataAnna/FT_bg_knm2golden_pseudoring_3.txt');

nbin(:,1)  = numel(FT_bg_knm2golden_pseudoring0(starfit:end,1));
hv(:,1)    = FT_bg_knm2golden_pseudoring0(starfit:end,1);
c(:,1)     = FT_bg_knm2golden_pseudoring0(starfit:end,2);
t(:,1)     = FT_bg_knm2golden_pseudoring0(starfit:end,3);

nbin(:,2)  = numel(FT_bg_knm2golden_pseudoring1(starfit:end,1));
hv(:,2)    = FT_bg_knm2golden_pseudoring1(starfit:end,1);
c(:,2)     = FT_bg_knm2golden_pseudoring1(starfit:end,2);
t(:,2)     = FT_bg_knm2golden_pseudoring1(starfit:end,3);

nbin(:,3)  = numel(FT_bg_knm2golden_pseudoring2(starfit:end,1));
hv(:,3)    = FT_bg_knm2golden_pseudoring2(starfit:end,1);
c(:,3)     = FT_bg_knm2golden_pseudoring2(starfit:end,2);
t(:,3)     = FT_bg_knm2golden_pseudoring2(starfit:end,3);

nbin(:,4)  = numel(FT_bg_knm2golden_pseudoring3(starfit:end,1));
hv(:,4)    = FT_bg_knm2golden_pseudoring3(starfit:end,1);
c(:,4)     = FT_bg_knm2golden_pseudoring3(starfit:end,2);
t(:,4)     = FT_bg_knm2golden_pseudoring3(starfit:end,3);

% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
psr  = 4;
Data = [hv(:,psr)/1e3 c(:,psr)./t(:,psr)*1e3 sqrt(c(:,psr))./t(:,psr)*1e3];

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
legend([hdata hfit],sprintf('Data - PSR %.0f',psr),sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(1),err(1),par(2),err(2)),'Location','SouthWest')
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
% Figure 1
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
hdata = semilogy(HVmin,parV(:,2)+2*errV(:,2),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('m+2\sigma (mcps/keV)','FontSize',14) 
grid on
PrettyFigureFormat;
publish_figurePDF(2,'FT_bg_knm2golden_2.pdf');

%% Case 3) Compare Pseudo-Rings
% Full Range
for psr=1:1:4
DataPSR = [hv(:,psr)/1e3 c(:,psr)./t(:,psr)*1e3 sqrt(c(:,psr))./t(:,psr)*1e3];
ParIni  = [mean(DataPSR(:,2)) 0]; 
Args    = {ParIni, DataPSR, '-c', 'min; minos'};
[ parPSR(psr,:), errPSR(psr,:), chi2minPSR(psr), errmatPSR(psr,:,:) ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline = %.3f ± %.3f\n',parPSR(psr,1),errPSR(psr,1));
fprintf('  slope    = %.3f ± %.3f\n',parPSR(psr,2),errPSR(psr,2));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));
fprintf('================================================================\n');
end
% Figure 3
figure(3)
subplot(3,1,2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
hdata = errorbar([1 2 3 4],parPSR(:,2),errPSR(:,2),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('PSR','FontSize',14) 
ylabel('slope (mcps/keV)','FontSize',14) 
PrettyFigureFormat;
grid on
xlim([0.5 4.5]);
subplot(3,1,1)
hdata = errorbar([1 2 3 4],parPSR(:,1),errPSR(:,1),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('PSR','FontSize',14) 
ylabel('baseline (mcps)','FontSize',14) 
PrettyFigureFormat;
grid on
xlim([0.5 4.5]);
subplot(3,1,3)
hdata = errorbar([1 2 3 4],parPSR(:,2)./parPSR(:,1).*sum(parPSR(:,1)),errPSR(:,2)./parPSR(:,1).*sum(parPSR(:,1)),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
xlabel('PSR','FontSize',14) 
ylabel('slope Uni.Eq (mcps/keV)','FontSize',14) 
PrettyFigureFormat;
grid on
xlim([0.5 4.5]);
publish_figurePDF(3,'FT_bg_knm2golden_3.pdf');

% Display Recommended Limit m+1sigma
SlopeLimit1sigma = table([1 2 3 4]',abs(errPSR(:,2)),norminv(0.68269,abs(parPSR(:,2)),errPSR(:,2)),'VariableNames',{'Pseudo-Ring','CurrentLimit (mcps/keV)','RecommendedLimit (mcps/keV)'});
SlopeLimit1sigma = varfun(@(var) round(var, 2), SlopeLimit1sigma)
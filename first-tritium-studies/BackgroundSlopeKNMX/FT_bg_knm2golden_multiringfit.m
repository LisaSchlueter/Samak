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
c(:,1)     = FT_bg_knm2golden_pseudoring0(starfit:end,2)*0.6811;
ce(:,1)    = sqrt(FT_bg_knm2golden_pseudoring0(starfit:end,2)*0.6811);
t(:,1)     = FT_bg_knm2golden_pseudoring0(starfit:end,3);

nbin(:,2)  = numel(FT_bg_knm2golden_pseudoring1(starfit:end,1));
hv(:,2)    = FT_bg_knm2golden_pseudoring1(starfit:end,1);
c(:,2)     = FT_bg_knm2golden_pseudoring1(starfit:end,2)*0.6811;
ce(:,2)    = sqrt(FT_bg_knm2golden_pseudoring1(starfit:end,2)*0.6811);
t(:,2)     = FT_bg_knm2golden_pseudoring1(starfit:end,3);

nbin(:,3)  = numel(FT_bg_knm2golden_pseudoring2(starfit:end,1));
hv(:,3)    = FT_bg_knm2golden_pseudoring2(starfit:end,1);
c(:,3)     = FT_bg_knm2golden_pseudoring2(starfit:end,2)*0.6811;
ce(:,3)    = sqrt(FT_bg_knm2golden_pseudoring2(starfit:end,2)*0.6811);
t(:,3)     = FT_bg_knm2golden_pseudoring2(starfit:end,3);

nbin(:,4)  = numel(FT_bg_knm2golden_pseudoring3(starfit:end,1));
hv(:,4)    = FT_bg_knm2golden_pseudoring3(starfit:end,1);
c(:,4)     = FT_bg_knm2golden_pseudoring3(starfit:end,2)*0.6811;
ce(:,4)    = sqrt(FT_bg_knm2golden_pseudoring3(starfit:end,2)*0.6811);
t(:,4)     = FT_bg_knm2golden_pseudoring3(starfit:end,3);

% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
Data = [reshape(hv,size(hv,1)*size(hv,2),1)/1e3  ...
        reshape(c./t*1e3,size(hv,1)*size(hv,2),1) ...
        reshape(ce./t*1e3,size(hv,1)*size(hv,2),1)];

% Linear Model
parnames = {'baseline1','slope1','baseline2','slope2','baseline3','slope3','baseline4','slope4'}; 
ParIni   = [mean(Data(:,2)) 0 mean(Data(:,2)) 0 mean(Data(:,2)) 0 mean(Data(:,2)) 0]; 

% Fit
Args = {ParIni, Data, '-c', 'min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2GaussMultiRing',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline1 = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope1    = %.3f ± %.3f - limit knm2 < %.1f mcps/keV \n',par(2),err(2),norminv(0.68269,abs(par(:,2)),err(:,2)));
fprintf('  baseline2 = %.3f ± %.3f\n',par(3),err(3));
fprintf('  slope2   = %.3f ± %.3f - limit knm2 < %.1f mcps/keV \n',par(4),err(4),norminv(0.68269,abs(par(:,4)),err(:,4)));
fprintf('  baseline3 = %.3f ± %.3f\n',par(5),err(5));
fprintf('  slope3    = %.3f ± %.3f - limit knm2 < %.1f mcps/keV \n',par(6),err(6),norminv(0.68269,abs(par(:,6)),err(:,6)));
fprintf('  baseline4 = %.3f ± %.3f\n',par(7),err(7));
fprintf('  slope4    = %.3f ± %.3f - limit knm2 < %.1f mcps/keV \n',par(8),err(8),norminv(0.68269,abs(par(:,8)),err(:,8)));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,sum(nbin)-numel(ParIni));
fprintf('================================================================\n');

% Figure 1
close all
figure(1)
set(gcf, 'Position',  [100, 100, 1400, 1200])

subplot(4,1,1)
hdata = errorbar(Data(1:26,1),Data(1:26,2),Data(1:26,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(1:26,1),Model(par(1:2),Data(1:26,1)),'LineWidth',3);
hold off
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(1),err(1),par(2),err(2)),'Location','EastOutside')
grid on

subplot(4,1,2)
hdata = errorbar(Data(27:52,1),Data(27:52,2),Data(27:52,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(27:52,1),Model(par(3:4),Data(27:52,1)),'LineWidth',3);
hold off
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(3),err(3),par(4),err(4)),'Location','EastOutside')
grid on

subplot(4,1,3)
hdata = errorbar(Data(53:78,1),Data(53:78,2),Data(53:78,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(53:78,1),Model(par(5:6),Data(53:78,1)),'LineWidth',3);
hold off
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(5),err(5),par(6),err(6)),'Location','EastOutside')
grid on

subplot(4,1,4)
hdata = errorbar(Data(79:104,1),Data(79:104,2),Data(79:104,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(79:104,1),Model(par(7:8),Data(79:104,1)),'LineWidth',3);
hold off;
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(7),err(7),par(8),err(8)),'Location','EastOutside')
grid on

publish_figurePDF(1,'./plots/FT_bg_knm2golden_multiring_1.pdf');


%% Figure 2
figure(2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
corplot(errmat);
PrettyFigureFormat;

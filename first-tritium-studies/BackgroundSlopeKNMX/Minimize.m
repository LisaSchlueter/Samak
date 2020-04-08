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

% Fit
Args = {ParIni, Data, '-c', 'min; minos , imp'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope    = %.3f ± %.3f\n',par(2),err(2));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));
fprintf('================================================================\n');

% Figure
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
legend([hdata hfit],'Data','Fit','Location','SouthWest')
grid on
publish_figurePDF(1,'FT_bg_knm2golden.pdf');

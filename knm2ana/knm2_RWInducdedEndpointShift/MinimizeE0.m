% Example of function minimization with fminuit
% Th. Lasserre, Feb 2013

% Code example with parabolic data
% The definition of chi2 (cost function to minimize) is in Chi2Poisson.m files
% The model definition is in model.m

addpath('./fminuit');
nbin = 2;

ParIni = 9.58; % Non-linear fit behaviour is dependent on initial guess => be careful!

b  = [1 2];
x  = [18573.66 18573.71]-18573.672;
ex = [0.008 0.014];
Data = [b' x' ex'];

parnames = 'mean'; 

% Fit
Args = {ParIni, Data, '-c', 'min; imp'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  mean = %.3f ± %.3f\n',par(1),err(1));
fprintf('  chi2min = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));

% Figure
figure(1)
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit = line([1,2],[Model(par,Data(:,1)) Model(par,Data(:,1))],'LineWidth',3);
hfiteu = line([1,2],[Model(par+err,Data(:,1)) Model(par+err,Data(:,1))],'LineWidth',1,'LineStyle','--');
hfited = line([1,2],[Model(par-err,Data(:,1)) Model(par-err,Data(:,1))],'LineWidth',1,'LineStyle','--');
hold off;
%xlabel('Data (1) 56160-56479         Data (2) 56560-56626','FontSize',24) 
ylabel('<E_0> - Cte (eV)','FontSize',24) 
            PrettyFigureFormat('FontSize',24); 
legend([hdata hfit],'KNM2 Data',sprintf('Fit pvalue = %.4f',chi2pvalue(chi2min,1)),'Location','SouthWest','FontSize',24)
grid off
xlim([0.9 2.1]);
xticks([1 2])
xticklabels({'56160-56479 ','56560-56626'})
legend('boxoff');

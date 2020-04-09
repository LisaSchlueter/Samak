% Example of function minimization with fminuit
% Th. Lasserre, Feb 2013

% Code example with parabolic data
% The definition of chi2 (cost function to minimize) is in Chi2Poisson.m files
% The model definition is in model.m

addpath('./fminuit');
nbin = 2;

ParIni = [9,1,1,1]; % Non-linear fit behaviour is dependent on initial guess => be careful!

b  = [1 2];
x  = [9.624 9.571];
ex = [0.012 0.031];
Data = [b' x' ex'];

parnames = 'mean'; 

% Fit
Args = {ParIni, Data, '-c', 'min; imp'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  mean = %.3f � %.3f\n',par(1),err(1));
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
xlabel('Data (1) using \tau (2) using \lambda','FontSize',14) 
ylabel('IBD prefactor (10^{-48} cm^2/MeV^2)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data','Fit','Location','SouthWest')
grid on
publish_figure(1,'fitibd.eps');

% % Example of function minimization with fminuit
% % Th. Lasserre, 2020
% 
Mainz.Experiment='Mainz';
Mainz.Reference='Eur.Phys.J. C40 (2005) 447-468';
Mainz.mBetaSquared=-0.6;
Mainz.Stat=2.2;
Mainz.Sys=2.1;
Mainz.Tot=sqrt(Mainz.Stat^2+Mainz.Sys^2);
Mainz.Marker='d';
Mainz.Color=rgb('CadetBlue');

Troitsk.Experiment='Troitsk';
Troitsk.Reference='Phys.Rev. D84 (2011) 112003';
Troitsk.mBetaSquared=-0.67;
Troitsk.Stat=1.89;
Troitsk.Sys=1.68;
Troitsk.Tot=sqrt(Troitsk.Stat^2+Troitsk.Sys^2);
Troitsk.Marker='s';
Troitsk.Color=rgb('IndianRed');


nbin = 2;

ParIni = 9.58; % Non-linear fit behaviour is dependent on initial guess => be careful!

b  = [1 2];
x  = [Mainz.mBetaSquared Troitsk.mBetaSquared];
ex = [Mainz.Tot Troitsk.Tot];
Data = [b' x' ex'];

parnames = 'mean'; 

% Fit
Args = {ParIni, Data, '-c', 'min; imp ; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  m^2 = %.3f ± %.3f\n',par(1),err(1));
fprintf('  chi2min = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));
fprintf('================================================================\n');

% Figure
figure(1)
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit = line([0,3],[Model(par,Data(:,1)) Model(par,Data(:,1))],'LineWidth',3);
hfiteu = line([0,3],[Model(par+err,Data(:,1)) Model(par+err,Data(:,1))],'LineWidth',1,'LineStyle','--');
hfited = line([0,3],[Model(par-err,Data(:,1)) Model(par-err,Data(:,1))],'LineWidth',1,'LineStyle','--');
hold off;
xlabel('(1) Mainz               (2) Troitsk','FontSize',14) 
ylabel('m_\nu^2','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'Data Mainz(1) Troitsk(2)',sprintf('Fit : m_\\nu^2 = %.2f \\pm %.2f eV^2',par,err),'Location','northoutside','FontSize',20);
grid on
xlim([0.5 2.5]);
export_fig(1,'MainzTroistkPull.pdf');

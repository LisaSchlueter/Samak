samakIS  = [0.795697049534547   0.174004647047996   0.026832295721868   0.003207320188824   0.000314360043257   0.000026209322691   0.000001906857616   0.000000123334562];
masterIS = [0.795622284380000   0.173996421603000   0.026831420781000   0.003207265664550   0.000314359541453   0.000026209669530   0.000001906910444   0.000000131449729];

samakIS   = samakIS(1:7);
masterIS  = masterIS(1:7);

ratioMasterSamak    = masterIS./samakIS;
absoluteDifference = (samakIS-masterIS);
relativeDifference = (samakIS-masterIS)./(samakIS+masterIS)/2*100;

% SuperImpose Master / Samak
fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Inelastic Scattering Probabilities - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G',...
    1.1e17,2.52,6.3);
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
r=bar((ratioMasterSamak-1));
set(gca, 'YScale', 'lin');
xlabel('Scattering','FontSize',18);
set(gca,'xticklabel',[0 1 2 3 4 5 6 7])
ylabel('Master / Samak -1','FontSize',18);
legend([r],'Master IS / Samak IS','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',18);

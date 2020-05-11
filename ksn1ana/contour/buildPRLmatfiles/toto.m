
file1      = importdata('coord_90eV_Twin_stat_95.mat');
file2      = importdata('coord_90eV_Twin_stat_95_freeM.mat');
file3      = importdata('coord_90eV_Twin_stat_95_fixM.mat');

figure(1)
h1=loglog(file1.sith4_X,file1.m4_Y,'LineWidth',2);
hold on
h2=loglog(file2.sith4_X,file2.m4_Y,'LineWidth',2);
h3=loglog(file3.sith4_X,file3.m4_Y,'LineWidth',2);
legend([h1 h2 h3],'m=0','m free','m^2=-1');
hold off
% Axis 
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('sin^2(2\theta_{ee})');
ylabel('\Deltam_{41}^2  (eV^2)');
PrettyFigureFormat

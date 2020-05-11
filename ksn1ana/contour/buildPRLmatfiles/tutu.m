close all
% 
% file1      = importdata('./contourmatfiles/coord_90eV_Real_stat.mat');
% file2      = importdata('./contourmatfiles/coord_90eV_Real_syst.mat');
% file3      = importdata('./contourmatfiles/coord_90eV_Real_syst_90_thierry.mat');
% 

file1      = importdata('./contourmatfiles/coord_95eV_Real_stat_90_thierry.mat');
file2      = importdata('./contourmatfiles/coord_95eV_Real_stat_95_thierry.mat');
file3      = importdata('./contourmatfiles/coord_95eV_Real_stat_99_thierry.mat');
file4      = importdata('./contourmatfiles/coord_95eV_Real_syst_90_thierry.mat');
file5      = importdata('./contourmatfiles/coord_95eV_Real_syst_95_thierry.mat');
file6      = importdata('./contourmatfiles/coord_95eV_Real_syst_99_thierry.mat');


figure(1)
h1=plot(file1.sith4_X,sqrt(file1.m4_Y),'LineWidth',2);
hold on
h2=plot(file2.sith4_X,sqrt(file2.m4_Y),'LineWidth',2);
h3=plot(file3.sith4_X,sqrt(file3.m4_Y),'LineWidth',2);
h4=plot(file4.sith4_X,sqrt(file4.m4_Y),'LineWidth',2);
h5=plot(file5.sith4_X,sqrt(file5.m4_Y),'LineWidth',2);
h6=plot(file6.sith4_X,sqrt(file6.m4_Y),'LineWidth',2);
legend([h1 h2 h3 h4 h5 h6],...
    'Thierry stat 90CL','Thierry stat 95CL','Thierry stat 99CL',...
    'Thierry stat+sys 90CL','Thierry stat+sys 95CL','Thierry stat+sys 99CL',...
    'Location','SouthWest');
hold off
% Axis 
grid on
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('sin^2(\theta_{ee})');
ylabel('m_{4}  (eV)');
axis([0.005 0.5 0.5 100])
PrettyFigureFormat

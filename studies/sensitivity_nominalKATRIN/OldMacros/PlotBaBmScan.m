%
% Plot KATRIN Sensitivity Scans in (Field Reduction / Ba) plane
% 
% Thierry, 11/10/2018
% 
%

% Define Range for MTD
range = 60;
% Load data
stat = load(sprintf('KATRIN_5y_BaBmScan_%0.fMTD_Stat.mat',range));
syst = load(sprintf('KATRIN_5y_BaBmScan_%0.fMTD_StatSyst.mat',range));

% Sensitivity Scan plot
figure(1)
x0=10;
y0=10;
width=1600;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])
Y=stat.Y./stat.Y(end:end);
subplot(1,2,1)
 [Cstat,hstat]=contourf(stat.X,Y,sqrt(stat.xxxSense)','Color','blue'); 
 hstat.LevelList = round(hstat.LevelList,2);  %rounds levels to 3rd decimal place
 hcl=clabel(Cstat,hstat,'FontSize',10,'Color','blue');
%surfc(stat.X,Y,sqrt(stat.xxxSense)');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'KATRIN m_\nu (90CL) - Preliminary')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('5 x 180 d - 60eV MTD - stat only');
PrettyFigureFormat
subplot(1,2,2)
 [Csyst,hsyst]=contourf(syst.X,Y,sqrt(syst.xxxSense)','Color','blue'); 
 hsyst.LevelList = round(hsyst.LevelList,2);  %rounds levels to 3rd decimal place
 hcl=clabel(Csyst,hsyst,'FontSize',10,'Color','blue');
%surfc(syst.X,Y,sqrt(syst.xxxSense)');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'KATRIN m_\nu (90CL) - Preliminary')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('5 x 180 d - 60eV MTD - stat + sys (samak)');
PrettyFigureFormat
publish_figurePDF(1,'Katrin5y60MTD_BaBmBsSensitivity_1.pdf')

% Systematics
figure(3)
x0=10;
y0=10;
width=1600;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])
Y=Y./Y(end:end);
subplot(1,2,1)
[Cstat,hstat]=contourf(stat.X,Y,sqrt(syst.xxxSense'-stat.xxxSense'),'Color','blue'); 
hstat.LevelList = round(hstat.LevelList,2);  %rounds levels to 3rd decimal place
hcl=clabel(Cstat,hstat,'FontSize',10,'Color','blue');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'KATRIN m_\nu (90CL)')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('5 x 180 d - 60eV MTD - stat only');
PrettyFigureFormat
subplot(1,2,2)
[Csyst,hsyst]=contourf(syst.X,Y,sqrt(syst.xxxSense'-stat.xxxSense')./sqrt(stat.xxxSense)','Color','blue'); 
hsyst.LevelList = round(hsyst.LevelList,2);  %rounds levels to 3rd decimal place
hcl=clabel(Csyst,hsyst,'FontSize',10,'Color','blue');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'KATRIN m_\nu (90CL)')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('5 x 180 d - 60eV MTD - stat + sys (samak)');
PrettyFigureFormat
publish_figurePDF(1,'Katrin5y60MTD_BaBmBsSensitivity_3.pdf')

% Background / Resolution Behaviohr
figure(2)
x0=10;
y0=10;
width=1600;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])
Y=Y./Y(end:end);
subplot(1,2,1)
[Cstat,hstat]=contourf(stat.X,Y,stat.B,'Color','blue'); 
hstat.LevelList = round(hstat.LevelList,2);  %rounds levels to 3rd decimal place
hcl=clabel(Cstat,hstat,'FontSize',10,'Color','blue');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'mcps - Preliminary')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('FPD Background Data+Simulation (mcps)');
PrettyFigureFormat
subplot(1,2,2)
[Csyst,hsyst]=contourf(syst.X,Y,stat.R,'Color','blue'); 
hsyst.LevelList = round(hsyst.LevelList,2);  %rounds levels to 3rd decimal place
hcl=clabel(Csyst,hsyst,'FontSize',10,'Color','blue');
colormap(flipud(pink));
hcb = colorbar; ylabel(hcb, 'eV')
xlabel('B_A (G)');
ylabel('B_m, Bs Reduction Factor');
title('MACE Resolution at 18575 eV');
PrettyFigureFormat
publish_figurePDF(2,'Katrin5y60MTD_BaBmBsSensitivity_2.pdf')

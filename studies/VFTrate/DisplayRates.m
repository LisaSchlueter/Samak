%
% Display Rates VFT
% Test Scan 19/05/2018
%

close all;

Z=ref_VFT(); Z.PlotTD;

% Compute Spectra Model
Z.ComputeTBDDS; Z.ComputeTBDIS;

% Data VFT1
%Count    = [9e4 3e4 1e4 1e3 1e2 0.5 0.15]'.*Z.TimeSec.*Z.qUfrac;
%Count    = [81289.4833 35939.8983 10660.8644 1207.338983 121.428571 0.5 0.179402]'.*Z.TimeSec.*Z.qUfrac;

% Data VFT2 - Trial
%Count    = [59345.7931      40516.1379      25929.6207      14910.3448      7620.65517      3079.57576        854.3125      333.586207       84.482759       51.103448              31       16.965517        8.159091        5.444444          4.0625        2.094595        1.392857        0.905263        0.609524        0.565217           0.256        0.404412        0.256579        0.353333        0.337748            0.26]';

% Data VFT2 - RUN 4257
%run40257 = importdata('data/spectrum_40257_outerRingsExcl.txt'); runcount=runcount+1;
run40257 = importdata('data/spectrum_40257.txt'); 
Count = flip(run40257(1:end-2,2));

% Data VFT2 - RUN 4258
%run40258 = importdata('data/spectrum_40258_outerRingsExcl.txt'); runcount=runcount+1;
%run40258 = importdata('data/spectrum_40258.txt'); runcount=runcount+1;
%Count = (run40258(1:end-2,2));

Count    = Count.*Z.TimeSec.*Z.qUfrac;
ErrCount = sqrt(Count);
Data = {Z.qU,Count,ErrCount};

%% Plot Data / Model
figure(2)
plt1= Plot(Z.qU-Z.Q,Z.TBDIS./(Z.TimeSec.*Z.qUfrac),Z.qU-Z.Q,Count./(Z.TimeSec.*Z.qUfrac));
%plt1= Plot(Z.qU-Z.Q,Z.TBDIS./(Z.TimeSec.*Z.qUfrac));
%plt1= Plot(Z.qU,Z.TBDIS./(Z.TimeSec.*Z.qUfrac));
plt1.LineWidth = 2;
plt1.LineStyle = {'-','--'};
plt1.Markers   = {'s','d'};
plt1.MarkerSpacing = [ 1 1];
pltstr         = sprintf('VFT Test Scan - 19/05/2018 - Preliminary');
plt1.Title     = pltstr; 
plt1.XLabel    = 'qU';  
plt1.YLabel    = 'cps';  
plt1.YScale    = 'log';  
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'Model (Samak)', 'Data'};
plt1.LegendLoc = 'NorthEast';
plttitle       = sprintf('Run Time %g sec',Z.TimeSec);
plt1.export('VFTtestScan1.png');

%% Plot Data / Model
figure(3)
plt1= Plot(Z.qU-Z.Q,Count./Z.TBDIS);
%plt1= Plot(Z.qU-Z.Q,Z.TBDIS./(Z.TimeSec.*Z.qUfrac));
%plt1= Plot(Z.qU,Z.TBDIS./(Z.TimeSec.*Z.qUfrac));
plt1.LineWidth = 2;
plt1.LineStyle = {'-'};
plt1.Markers   = {'s'};
plt1.MarkerSpacing = [ 1];
pltstr         = sprintf('VFT Test Scan');
plt1.Title     = pltstr; 
plt1.XLabel    = 'qU';  
plt1.YLabel    = 'Data/Model';  
plt1.YScale    = 'log';  plt1.YLim=[0,1];
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'Samak Model'};
plt1.LegendLoc = 'SouthWest';
plttitle       = sprintf('Run Time %g sec',Z.TimeSec);
plt1.export('VFTtestScan2.png');

%% Plot With Errorbars
figure(10)
subplot(2,1,1)
%h1=errorbar(Z.qU-Z.Q,Count./(Z.TimeSec.*Z.qUfrac),ErrCount./(Z.TimeSec.*Z.qUfrac),'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
h1=errorbar(Z.qU-Z.Q,Count,ErrCount,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
hold on
h2=plot(Z.qU-Z.Q,Z.TBDIS,'LineStyle','--','LineWidth',1);
hold off
set(gca, 'YScale', 'log');
grid on
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('cps','FontSize',14);
title('VFT - 19/05/2018 - Preliminary')
set(gca,'FontSize',12);
l1 = sprintf('Data');
l2 = sprintf('Samak Model (no fit)');
a = legend([h1 h2],l1,l2);
PrettyFigureFormat
subplot(2,1,2)
errorbar(Z.qU-Z.Q,Count./Z.TBDIS,ErrCount./Z.TBDIS,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('ratio','FontSize',14);
grid on
PrettyFigureFormat
set(gcf, 'Position', [100, 100, 1000, 500])
saveas(gcf,'vft-dataVmodel.png')
publish_figure(10,'vft-dataVmodel.eps');

%% Plot With Errorbars
figure(10)
subplot(2,1,1)
%h1=errorbar(Z.qU-Z.Q,Count./(Z.TimeSec.*Z.qUfrac),ErrCount./(Z.TimeSec.*Z.qUfrac),'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
h1=errorbar(Z.qU-Z.Q,Count,ErrCount,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
hold on
h2=plot(Z.qU-Z.Q,Z.TBDIS,'LineStyle','--','LineWidth',1);
hold off
set(gca, 'YScale', 'log');
grid on
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('cps','FontSize',14);
title('VFT - 19/05/2018 - Preliminary')
set(gca,'FontSize',12);
l1 = sprintf('Data');
l2 = sprintf('Samak Model (no fit)');
a = legend([h1 h2],l1,l2);
PrettyFigureFormat
subplot(2,1,2)
h1=errorbar(Z.qU-Z.Q,Count,ErrCount,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
hold on
h2=plot(Z.qU-Z.Q,Z.TBDIS,'LineStyle','--','LineWidth',1);
hold off
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('counts','FontSize',14);
grid on
PrettyFigureFormat
set(gcf, 'Position', [100, 100, 1000, 500])
saveas(gcf,'vft-dataVmodel.png')
publish_figure(10,'vft-dataVmodel.eps');


%%
for i=1:1:Z.nqU
    fprintf(2,'qU=%f V - Rate=%f cps\n',Z.qU(i),Z.TBDIS(i)./(Z.TimeSec.*Z.qUfrac(i)))
end
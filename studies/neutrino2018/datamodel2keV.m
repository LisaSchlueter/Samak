addpath(genpath('../../../Samak2.0'));
close all

run = 40263;

%% Settings

% Lower  Bin 9  = -200 V
FitLB = 1;
% Higher Bin 26 = + 40 V
FitHB = 26;

chi2 = 'Stat'; %CM, CMFrac, Stat
RunNr = 40263;
exclDataStart=FitLB; % Start Fit at bin i (i=1 -> qU(min)=-1.7keV, i=9 ->qU(min)=200eV)

%% Model
M = ref_runsummaries(RunNr,'ISCS','Theory','recomputeRF','OFF'); 
M.ComputeTBDDS; M.ComputeTBDIS;
CountSim    = M.TBDIS; CountSimErr = sqrt(M.TBDIS); 
Sim  = [M.qU,CountSim,CountSimErr];

%% Read Data
D = importdata([num2str(RunNr),'.mat']);
switch RunNr
    case {40259,40260,40263, 40264, 40265,40266}
        Data = [D.qU, D.TBDIS, sqrt(D.TBDIS)];
    case {40257,40258}
        Data = [D.qU(3:end), D.TBDIS(3:end), sqrt(D.TBDIS(3:end))];
end
    
%% Covariance Matrix
cm = importdata('CovMat_Run40257_0.1rhod_0.02BField_0.02ISX_0.03FSD_0.023TASR_TC.mat');
covmatfrac = cm.CovMatFracCombi;
CM_Stat = diag(Data(:,2));
FitCM = Data(:,2).*covmatfrac.*Data(:,2)' + CM_Stat;

%% Plot With Errorbars

%% Plot With Errorbars
figure(11)
set(gcf, 'Position', [100, 100, 1400, 800])
h1=errorbar(M.qU-M.Q,Data(:,2)./M.qUfrac./M.TimeSec,Data(:,3)./M.qUfrac./M.TimeSec,'ko','MarkerSize',3,'MarkerEdgeColor' , [0 0 0]);
h1.MarkerSize = 3; h1.CapSize = 0;h1.LineStyle= 'none';h1.LineWidth= 2;legend hide
hold on
hfit1 = boundedline(M.qU-M.Q,Sim(:,2)./M.qUfrac./M.TimeSec,diag(sqrt(FitCM))./M.qUfrac./M.TimeSec,'alpha','cmap',rgb('CadetBlue'));
%h2=plot(M.qU-M.Q,Sim(:,2),'LineStyle','--','LineWidth',3);
%h1=errorbar(M.qU-M.Q,Data(:,2)./M.qUfrac./M.TimeSec,Data(:,3)./M.qUfrac./M.TimeSec,'ko','MarkerSize',3,'MarkerEdgeColor' , [0 0 0]);
%h1.CapSize = 0; h1.MarkerEdgeColor = [0 0 0];
[ll la]= boundedline(M.qU-M.Q,Sim(:,2)./M.qUfrac./M.TimeSec,sqrt(diag(FitCM))./M.qUfrac./M.TimeSec,'alpha','cmap',rgb('CadetBlue'));
ll.LineStyle='--';
hold off
set(gca, 'YScale', 'lin');
%grid on
xlim([-1600 50])
xlabel('retarding potential -18575 (Volt)','FontSize',16);
ylabel('counts per second','FontSize',16);
set(gca,'FontSize',12);
l1 = sprintf('Data (Run %g, %g seconds)',run,M.TimeSec);
l2 = sprintf('Model and uncorrelated error band (no fit)');
a = legend([h1 la],l1,l2,'Location','NorthEast'); a.String=a.String(1:2);
legend('boxoff'); a = legend([h1 la],l1,l2,'Location','NorthEast');
axis([min(M.qU-M.Q) max(M.qU-M.Q) -1000 max(Data(:,2)./M.qUfrac./M.TimeSec)*1.1 ]);
title('KATRIN Tritium Commissioning Run - 19 May 2018');
PrettyFigureFormat

% zoomPlot to highlight a portion of the major plot 
[p,z] = zoomPlotError(M.qU-M.Q,Data(:,2)./M.qUfrac./M.TimeSec,Data(:,3)./M.qUfrac./M.TimeSec,[-220 40],[0.4 0.45 0.45 0.3],[3 4]); 
hold on 
set(gca, 'YScale', 'log');
[ll la]=boundedline(M.qU-M.Q,Sim(:,2)./M.qUfrac./M.TimeSec,sqrt(diag(FitCM))./M.qUfrac./M.TimeSec,'alpha','cmap',rgb('CadetBlue'));
ll.LineStyle='-';ll.LineWidth=1;

hold off
a = legend([h1 la],l1,l2,'Location','NorthEast');
PrettyFigureFormat
saveas(gcf,'figures/vft-dataVmodel_stat.png')
publish_figure(11,'figures/vft-dataVmodel_stat2.eps');
%title('KATRIN - First Tritium Test Run - 19/05/2018 - Preliminary')

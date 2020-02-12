addpath(genpath('../../../Samak2.0')); close all;
%% Settings
chi2 = 'CM'; %CM, CMFrac, Stat, CMShape
RunNr = 40263;
exclDataStart=9; % Start Fit at bin i (i=1 -> qU(min)=-1.7keV, i=9 ->qU(min)=200eV)
%% Read Data
D = importdata([num2str(RunNr),'.mat']);
switch RunNr
    case {40259,40260,40263, 40264, 40265,40266}
        ErrCount = sqrt(D.TBDIS(1:end));
        if RunNr == 40263
            %Count(4) = 99999;
            %ErrCount(4) = 99999;
        end
        Data = [D.qU, D.TBDIS, ErrCount];
    case {40257,40258}
        Data = [D.qU(3:end), D.TBDIS(3:end), sqrt(D.TBDIS(3:end))];
end
%% Model
M = ref_runsummaries(RunNr,'ISCS','Theory','recomputeRF','OFF');
M.ComputeTBDDS; M.ComputeTBDIS;

%% Covariance Matrix
%cm = importdata('CovMat_Run40257_0.1rhod_0.02BField_0.02ISX_0.03FSD_0.023TASR_TC.mat');
cm = importdata('CovMat_Run40257_1000-Trials_-RF_EL-RF_BF-RF_RX-FSD-TASR-TCoff_RAD-TCoff_OTHER.mat');
covmatfrac = cm.CovMatFracCombi;
CM_Stat = diag(Data(:,3).^2);
CMFrac_Stat = diag(1./Data(:,2));
FitCM = Data(:,2).*covmatfrac.*Data(:,2)' + CM_Stat;
FitCMFrac = covmatfrac + CMFrac_Stat; 
FitCMShape =  Data(:,2).*cm.obj.CovMatFracShape.*Data(:,2)' + CM_Stat;
FitCMNorm =  Data(:,2).*cm.obj.CovMatFracNorm.*Data(:,2)' + CM_Stat;
%FitCMsum =  Data(:,2).*(cm.obj.CovMatFracNorm+cm.obj.CovMatFracShape).*Data(:,2)' + CM_Stat;
%% Fit
tic;
switch chi2
    case 'CM'
%         F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
%             'chi2name','chi2CM', 'COVMAT', FitCM,...
%             'pulls',[4 inf inf inf],...
%             'exclDataStart',exclDataStart); % excludes start-1 points
                       F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
                    'chi2name','chi2CMShape', 'COVMATShape', FitCM,...
                    'pulls',[4 inf inf inf],...
                    'exclDataStart',exclDataStart); % excludes start-1 points
    case 'CMFrac'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
            'chi2name','chi2CMFrac', 'COVMATFrac', FitCMFrac,...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart);
    case 'Stat'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
            'chi2name','chi2',...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart);
            FitCM = CM_Stat;
    case 'CMShape'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...         
            'chi2name','chi2CM', 'COVMAT', FitCMShape,...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart);
        FitCM = abs(FitCMShape);
end
 toc;

par = F.RESULTS{1};
err = F.RESULTS{2};
chi2min = F.RESULTS{3};
dof = numel(M.qU(exclDataStart:end))-4;

%% Plot Spectrum + Fit
close 
fig1 = figure(5);
set(fig1, 'Units', 'normalized', 'Position', [0.9, 0.9, 0.8, 1.2]);
subplot(2,1,1);
[lfit pfit] = boundedline(M.qU-M.Q_i,M.TBDIS./M.qUfrac./M.TimeSec,diag(sqrt(FitCM)./M.qUfrac./M.TimeSec),'alpha','cmap',rgb('CadetBlue'));
hold on;
hdata = errorbar(Data(:,1)-M.Q_i,Data(:,2)./M.qUfrac./M.TimeSec,Data(:,3)./M.qUfrac./M.TimeSec,'ko','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hdata.CapSize = 0;
hold off;
xlim([min(M.qU(exclDataStart:end))-M.Q_i max(M.qU)-M.Q_i]);
myfitleg = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.2f \t\\pm %.2f eV \n E0=%.2f \t\\pm %.2f eV \n B=%.2f \t\\pm %.2f mcps\n N=%.2f \t\\pm %.2f',...
    chi2min,dof,((M.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
    (M.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,par(4)+1,err(4));
%leg1 = legend([hdata lfit],'Data',myfitleg,'Location','SouthWest') ;
legdata = sprintf('Data (Run %u, %u seconds)',RunNr,M.TimeSec);
% if strcmp(chi2,'CM')
%     chi2leg = 'with Covariance Matrix'
% end
%legfit = sprintf('Fit (%s)',chi2leg);
legend([hdata lfit],legdata, myfitleg,'Location','southwest');
xlabel('retarding potential - 18575 (V)');
ylabel('counts per second');
set(gca,'yscale','log');
%title(sprintf('Run %u Data - Samak Fit (%s) \n(18575 eV - qU_{min}) = %.0f eV',RunNr,chi2,M.Q_i-Data(exclDataStart,1)));
title(sprintf('KATRIN Tritium Commissioning Run - 19 May 2018'));
legend('boxoff')
PrettyFigureFormat;
set(gca,'YMinorTick','off');
set(gca,'TickLength',[0.01 0.01]);
ylim([0.75*min(M.TBDIS./M.qUfrac./M.TimeSec) 1.25*M.TBDIS(exclDataStart)./M.qUfrac(exclDataStart)./M.TimeSec]);

subplot(2,1,2);
[lstat pstat]  = boundedline(Data(:,1)-M.Q_i,zeros(26,1),sqrt(diag(CM_Stat))./sqrt(diag(FitCM)),'alpha','transparency', 0.4,'cmap',rgb('SteelBlue'));%[0 0 0]);%0.8*[0 0 0]);
hold on;
[lsys,psys] = boundedline(Data(:,1)-M.Q_i,zeros(26,1),ones(26,1),'alpha','transparency',0.4,'cmap',rgb('DarkCyan'));
lstat.LineStyle= '--'; lsys.LineStyle= '--';
plot(Data(:,1)-M.Q_i,(Data(:,2)-M.TBDIS)./sqrt(diag(FitCM)),'o','MarkerEdgeColor','k');
leg1 = legend([pstat psys],'stat','stat+sys','Location','NorthWest'); %hsyst
legend('boxoff')
hold off;
xlabel('retarding potential - 18575 (V)');
ylabel('norm. residuals');
PrettyFigureFormat;
xlim([min(M.qU(exclDataStart:end))-M.Q_i max(M.qU)-M.Q_i]);
ylim([1.2*min((Data(:,2)-M.TBDIS)./sqrt(diag(FitCM)))  2]);%2*max((Data(:,2)-M.TBDIS)./sqrt(diag(FitCM)))])
set(gca,'YMinorTick','off');
set(gca,'XMinorTick','off');
%set(gca,'YTick',[-1 1]);
set(gca,'TickLength',[0.01 0.01]);
save_name = sprintf('./plots/FitRun%u_chi2-%s-excl%u.eps',RunNr,chi2,exclDataStart);
publish_figure(fig1,save_name);

%
% FitVFT1
% 19/05/18
% T. Lasserre
%

clear all

% Model
A = ref_VFT;
A.ComputeTBDDS();
A.ComputeTBDIS();

% Data VFT1
%Count    = [81289.4833 35939.8983 10660.8644 1207.338983 121.428571 0.5 0.179402]'.*A.TimeSec.*A.qUfrac;
% Data VFT2
Count    = [59345.7931      40516.1379      25929.6207      14910.3448      7620.65517      3079.57576        854.3125      333.586207       84.482759       51.103448              31       16.965517        8.159091        5.444444          4.0625        2.094595        1.392857        0.905263        0.609524        0.565217           0.256        0.404412        0.256579        0.353333        0.337748            0.26]';
Count    = Count.*A.TimeSec.*A.qUfrac;

ErrCount = sqrt(Count)*1;
Data = {A.qU,Count,ErrCount};

% Covariance Matrix
myEffects = struct(...
    'RF_EL','ON',...   % Response Function(RF) EnergyLoss
    'RF_BF','ON',...   % RF B-Fields
    'RF_RX','ON',...   % RF Column Density, Cross Section
    'FSD','ON',...     % Activity Fluctuations
    'TCoff','ON',...   % Theoretical Corrections
    'TASR','ON');      % FSD

% C = CovarianceMatrix('StudyObject',A,...
%     'nTrials',1000,...
%     'RunFlag','OFF', 'nRuns',1,...
%     'SysEffect',myEffects,...
%     'MACE_Ba_T_RelErr',1e-2,...
%     'MACE_Bmax_T_RelErr',1e-2,...
%     'WGTS_B_T_RelErr',1e-2,...
%     'WGTS_CD_MolPerCm2_RelErr',0.1,...
%     'RecomputeFlag','OFF');
% % C.ComputeCM;
% a=load('MyCovMat.mat')
% CovMat = a.covmat;
% CovMat = nearestSPD(CovMat)+diag(Count)

% Fit
%F = FITC('SO',A,'DATA',Data,'COVMAT',CovMat,'fitter','matlab','chi2name','chi2');
F = FITC('SO',A,'DATA',Data,'fitter','minuit','chi2name','chi2');
F.RESULTS{1}
F.RESULTS{2}
F.RESULTS{3}

%% Plot Data / Fit
% Plot With Errorbars
figure(10)
subplot(2,1,1)
h1=errorbar(A.qU-A.Q,Count./(A.TimeSec.*A.qUfrac),ErrCount./(A.TimeSec.*A.qUfrac),'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
hold on
h2=plot(A.qU-A.Q,A.TBDIS./(A.TimeSec.*A.qUfrac),'LineStyle','--','LineWidth',1);
hold off
set(gca, 'YScale', 'log');
grid on
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('cps','FontSize',14);
title('VFT - 19/05/2018 - Preliminary')
set(gca,'FontSize',12);
l1 = sprintf('Data');
l2 = sprintf('Samak Model (fit)');
a = legend([h1 h2],l1,l2);
PrettyFigureFormat
subplot(2,1,2)
errorbar(A.qU-A.Q,Count./A.TBDIS,ErrCount./A.TBDIS,'ks','MarkerSize',3,'MarkerFaceColor',.8*[1 1 1],'LineWidth',3);
xlabel('qU-E_0 (V)','FontSize',14);
ylabel('ratio','FontSize',14);
grid on
PrettyFigureFormat
set(gcf, 'Position', [100, 100, 1000, 500])
%saveas(gcf,'vft-dataVmodel.png')
publish_figure(10,'vft-dataVmodel.eps');

% 

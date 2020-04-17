% 
% First Tritium Background Slope Constraint
% Th. Lasserre, April 2020
%

%% Read KNM1 & KNM2

% Read KNM1 data
range     = 5;
RunAnaArgKNM1 = {'RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
    'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
    'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
    'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
    'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
    'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
    'minuitOpt','min ; minos',...         % technical fitting options (minuit)
    'FSDFlag','Sibille0p5eV',...          % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1',...
    'Twin_SameCDFlag','ON',...
    'Twin_SameIsotopFlag','ON',...
    'RingMerge','Full'...
    };%
NonPoissonScaleFactorKNM1 = 1.064;
KNM1                      = MultiRunAnalysis(RunAnaArgKNM1{:});
KNM1.exclDataStart        = KNM1.GetexclDataStart(range);
KNM1PSR                   = RingAnalysis('RunAnaObj',KNM1,'RingList',1:4);
KNM1PSR1                  = KNM1PSR.MultiObj(1);
KNM1PSR2                  = KNM1PSR.MultiObj(2);
KNM1PSR3                  = KNM1PSR.MultiObj(3);
KNM1PSR4                  = KNM1PSR.MultiObj(4);

% Read KNM2 data
range         = 4;
RunAnaArgKNM2 = {...
    'RunList','KNM2_Prompt',...        % all KNM2 golden runs
    'DataType','Real',...
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'ROIFlag','Default',...
    'RingMerge','Full'...
    };%
NonPoissonScaleFactorKNM2 = 1.112;
KNM2                      = MultiRunAnalysis(RunAnaArgKNM2{:});
KNM2.exclDataStart        = KNM2.GetexclDataStart(range);
KNM2PSR                   = RingAnalysis('RunAnaObj',KNM2,'RingList',1:4);
KNM2PSR1                  = KNM2PSR.MultiObj(1);
KNM2PSR2                  = KNM2PSR.MultiObj(2);
KNM2PSR3                  = KNM2PSR.MultiObj(3);
KNM2PSR4                  = KNM2PSR.MultiObj(4);

%% Extract KNM1 rates
knm1psr1_nb   = numel(KNM1PSR1.RunData.qU(KNM1.exclDataStart:end));
knm1psr1_hv   = KNM1PSR1.RunData.qU(KNM1.exclDataStart:end);
knm1psr1_r    = KNM1PSR1.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1PSR1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR1.RunData.TimeSec;
knm1psr1_re   = NonPoissonScaleFactorKNM1*KNM1PSR1.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1PSR1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR1.RunData.TimeSec;

knm1psr2_nb   = numel(KNM1PSR2.RunData.qU(KNM1.exclDataStart:end));
knm1psr2_hv   = KNM1PSR2.RunData.qU(KNM1.exclDataStart:end);
knm1psr2_r    = KNM1PSR2.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1PSR2.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR2.RunData.TimeSec;
knm1psr2_re   = NonPoissonScaleFactorKNM1*KNM1PSR2.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1PSR2.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR2.RunData.TimeSec;

knm1psr3_nb   = numel(KNM1PSR3.RunData.qU(KNM1.exclDataStart:end));
knm1psr3_hv   = KNM1PSR3.RunData.qU(KNM1.exclDataStart:end);
knm1psr3_r    = KNM1PSR3.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1PSR3.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR3.RunData.TimeSec;
knm1psr3_re   = NonPoissonScaleFactorKNM1*KNM1PSR3.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1PSR3.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR3.RunData.TimeSec;

knm1psr4_nb   = numel(KNM1PSR4.RunData.qU(KNM1.exclDataStart:end));
knm1psr4_hv   = KNM1PSR4.RunData.qU(KNM1.exclDataStart:end);
knm1psr4_r    = KNM1PSR4.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1PSR4.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR4.RunData.TimeSec;
knm1psr4_re   = NonPoissonScaleFactorKNM1*KNM1PSR4.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1PSR4.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1PSR4.RunData.TimeSec;

%% Extract KNM2 rates
knm2psr1_nb   = numel(KNM2PSR1.RunData.qU(KNM2PSR1.exclDataStart:end));
knm2psr1_hv   = KNM2PSR1.RunData.qU(KNM2PSR1.exclDataStart:end);
knm2psr1_r    = KNM2PSR1.RunData.TBDIS(KNM2PSR1.exclDataStart:end)./KNM2PSR1.RunData.qUfrac(KNM2PSR1.exclDataStart:end)/KNM2PSR1.RunData.TimeSec;
knm2psr1_re   = NonPoissonScaleFactorKNM2*KNM2PSR1.RunData.TBDISE(KNM2PSR1.exclDataStart:end)./KNM2PSR1.RunData.qUfrac(KNM2PSR1.exclDataStart:end)/KNM2PSR1.RunData.TimeSec;

knm2psr2_nb   = numel(KNM2PSR2.RunData.qU(KNM2PSR2.exclDataStart:end));
knm2psr2_hv   = KNM2PSR2.RunData.qU(KNM2PSR2.exclDataStart:end);
knm2psr2_r    = KNM2PSR2.RunData.TBDIS(KNM2PSR2.exclDataStart:end)./KNM2PSR2.RunData.qUfrac(KNM2PSR2.exclDataStart:end)/KNM2PSR2.RunData.TimeSec;
knm2psr2_re   = NonPoissonScaleFactorKNM2*KNM2PSR2.RunData.TBDISE(KNM2PSR2.exclDataStart:end)./KNM2PSR2.RunData.qUfrac(KNM2PSR2.exclDataStart:end)/KNM2PSR2.RunData.TimeSec;

knm2psr3_nb   = numel(KNM2PSR3.RunData.qU(KNM2PSR3.exclDataStart:end));
knm2psr3_hv   = KNM2PSR3.RunData.qU(KNM2PSR3.exclDataStart:end);
knm2psr3_r    = KNM2PSR3.RunData.TBDIS(KNM2PSR3.exclDataStart:end)./KNM2PSR3.RunData.qUfrac(KNM2PSR3.exclDataStart:end)/KNM2PSR3.RunData.TimeSec;
knm2psr3_re   = NonPoissonScaleFactorKNM2*KNM2PSR3.RunData.TBDISE(KNM2PSR3.exclDataStart:end)./KNM2PSR3.RunData.qUfrac(KNM2PSR3.exclDataStart:end)/KNM2PSR3.RunData.TimeSec;

knm2psr4_nb   = numel(KNM2PSR4.RunData.qU(KNM2PSR4.exclDataStart:end));
knm2psr4_hv   = KNM2PSR4.RunData.qU(KNM2PSR4.exclDataStart:end);
knm2psr4_r    = KNM2PSR4.RunData.TBDIS(KNM2PSR4.exclDataStart:end)./KNM2PSR4.RunData.qUfrac(KNM2PSR4.exclDataStart:end)/KNM2PSR4.RunData.TimeSec;
knm2psr4_re   = NonPoissonScaleFactorKNM2*KNM2PSR4.RunData.TBDISE(KNM2PSR4.exclDataStart:end)./KNM2PSR4.RunData.qUfrac(KNM2PSR4.exclDataStart:end)/KNM2PSR4.RunData.TimeSec;

%%

% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
k2i    = [1 ; 7];
hv     = [knm1psr1_hv ; knm2psr1_hv(k2i) ; knm1psr2_hv ; knm2psr2_hv(k2i) ; knm1psr3_hv ; knm2psr3_hv(k2i) ; knm1psr4_hv ; knm2psr4_hv(k2i)];
c      = [knm1psr1_r ; knm2psr1_r(k2i) ; knm1psr2_r ; knm2psr2_r(k2i) ; knm1psr3_r ; knm2psr3_r(k2i) ; knm1psr4_r ; knm2psr4_r(k2i)];
ce     = [knm1psr1_re ; knm2psr1_re(k2i) ; knm1psr2_re ; knm2psr2_re(k2i) ; knm1psr3_re ; knm2psr3_re(k2i) ; knm1psr4_re ; knm2psr4_re(k2i)];

Data = [hv/1e3 c*1e3 ce*1e3];

%%
% Linear Model
parnames = {'bknm1psr1','sknm1psr1','bknm1psr2','sknm1psr2','bknm1psr3','sknm1psr3','bknm1psr4','sknm1psr4','bknm2psr1','sknm2psr1','bknm2psr2','sknm2psr2','bknm2psr3','sknm2psr3','bknm2psr4','sknm2psr4'}; 
ParIni   = [mean(Data(1:7,2)) 0 mean(Data(8:9,2)) 0 mean(Data(10:16,2)) 0 mean(Data(17:18,2)) 0 mean(Data(19:25,2)) 0 mean(Data(26:27,2)) 0 mean(Data(28:34,2)) 0 mean(Data(35:36,2)) 0]; 

% Fit1
Args = {ParIni, Data, '-c', 'min; minos; imp'};
[ par, err, chi2min, errmat ] = fminuit('Chi2GaussKNM12MultiRing',Args{:});
% Fit2
 Args = {ParIni, Data, '-c', 'fix 1 3 5 7 ; fix 9 11 13 15; min; minos; imp'};
 [ par, err, chi2min, errmat ] = fminuit('Chi2GaussKNM12MultiRing',Args{:});

fprintf('=========================== Fit results ========================\n');

fprintf('  baseline knm1 psr1 = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope knm1 psr1    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(2),err(2),norminv(0.68269,abs(par(:,2)),err(:,2)));
fprintf('  baseline knm2 psr1 = %.3f ± %.3f\n',par(3),err(3));
fprintf('  slope knm2 psr1    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(4),err(4),norminv(0.68269,abs(par(:,4)),err(:,4)));

fprintf('  baseline knm1 psr2 = %.3f ± %.3f\n',par(5),err(5));
fprintf('  slope knm1 psr2    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(6),err(6),norminv(0.68269,abs(par(:,6)),err(:,6)));
fprintf('  baseline knm2 psr2 = %.3f ± %.3f\n',par(7),err(7));
fprintf('  slope knm2 psr2    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(8),err(8),norminv(0.68269,abs(par(:,8)),err(:,8)));

fprintf('  baseline knm1 psr3 = %.3f ± %.3f\n',par(9),err(9));
fprintf('  slope knm1 psr3    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(10),err(10),norminv(0.68269,abs(par(:,10)),err(:,10)));
fprintf('  baseline knm2 psr3 = %.3f ± %.3f\n',par(11),err(11));
fprintf('  slope knm2 psr3    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(12),err(12),norminv(0.68269,abs(par(:,12)),err(:,12)));

fprintf('  baseline knm1 psr4 = %.3f ± %.3f\n',par(13),err(13));
fprintf('  slope knm1 psr4    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(14),err(14),norminv(0.68269,abs(par(:,14)),err(:,14)));
fprintf('  baseline knm2 psr4 = %.3f ± %.3f\n',par(15),err(15));
fprintf('  slope knm2 psr4    = %.3f ± %.3f - limit < %.1f mcps/keV \n',par(16),err(16),norminv(0.68269,abs(par(:,16)),err(:,16)));

fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,numel(Data(:,1))-numel(ParIni));
fprintf('================================================================\n');

%ParIni   = [mean(Data(1:7,2)) 0 mean(Data(8:9,2)) 0 mean(Data(10:16,2)) 0 mean(Data(17:18,2)) 0 mean(Data(19:25,2)) 0 mean(Data(26:27,2)) 0 mean(Data(28:34,2)) 0 mean(Data(35:36,2)) 0]; 

%% Figure 1
myMainTitle = sprintf('KATRIN KNM1+KNM2 PSR-wise Background Slope Constraint - PSR-Uncorrelated - KNM12-Correlated');
maintitle   = myMainTitle;
fig1      = figure('Name',sprintf('KATRIN KNM1+KNM2 PSR-wise Background Slope Constraint - PSR-Corrleated - KNM12-Correlated'),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1200]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

subplot(4,2,1)
hdata = errorbar(Data(1:7,1),Data(1:7,2),Data(1:7,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on;
hfit    = plot(Data(1:7,1),Model(par(1:2),Data(1:7,1)),'LineWidth',3,'Color',rgb('IndianRed'));
hold off
xlabel('KNM1 Retarding Potential PSR1 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(2),err(2)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,2)
hdata = errorbar(Data(8:9,1),Data(8:9,2),Data(8:9,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2);
hold on;
hfit    = plot(Data(8:9,1),Model(par(3:4),Data(8:9,1)),'LineWidth',3,'Color',rgb('DarkBlue'));
hold off
xlabel('KNM2 Retarding Potential PSR1 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(4),err(4)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,3)
hdata = errorbar(Data(10:16,1),Data(10:16,2),Data(10:16,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on;
hfit    = plot(Data(10:16,1),Model(par(5:6),Data(10:16,1)),'LineWidth',3,'Color',rgb('IndianRed'));
hold off
xlabel('KNM1 Retarding Potential PSR2 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(6),err(6)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,4)
hdata = errorbar(Data(17:18,1),Data(17:18,2),Data(17:18,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2);
hold on;
hfit    = plot(Data(17:18,1),Model(par(7:8),Data(17:18,1)),'LineWidth',3,'Color',rgb('DarkBlue'));
hold off
xlabel('KNM2 Retarding Potential PSR2 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(8),err(8)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,5)
hdata = errorbar(Data(19:25,1),Data(19:25,2),Data(19:25,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on;
hfit    = plot(Data(19:25,1),Model(par(9:10),Data(19:25,1)),'LineWidth',3,'Color',rgb('IndianRed'));
hold off
xlabel('KNM1 Retarding Potential PSR3 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(10),err(10)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,6)
hdata = errorbar(Data(26:27,1),Data(26:27,2),Data(26:27,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2);
hold on;
hfit    = plot(Data(26:27,1),Model(par(11:12),Data(26:27,1)),'LineWidth',3,'Color',rgb('DarkBlue'));
hold off
xlabel('KNM2 Retarding Potential PSR3 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(12),err(12)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,7)
hdata = errorbar(Data(28:34,1),Data(28:34,2),Data(28:34,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on;
hfit    = plot(Data(28:34,1),Model(par(13:14),Data(28:34,1)),'LineWidth',3,'Color',rgb('IndianRed'));
hold off
xlabel('KNM1 Retarding Potential PSR4 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(14),err(14)),'Location','SouthEast','FontSize',8);
grid on

subplot(4,2,8)
hdata = errorbar(Data(35:36,1),Data(35:36,2),Data(35:36,3),'ks','MarkerSize',8,'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2);
hold on;
hfit    = plot(Data(35:36,1),Model(par(15:16),Data(35:36,1)),'LineWidth',3,'Color',rgb('DarkBlue'));
hold off
xlabel('KNM2 Retarding Potential PSR4 (keV)','FontSize',12) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hfit],sprintf('slope=%0.3f \\pm %0.3f mcps/keV',par(16),err(16)),'Location','SouthEast','FontSize',8);
grid on

publish_figurePDF(1,'./BKG_KNM1and2PSR.pdf');

return;


%% Figure 2
figure(2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
corplot(errmat);
PrettyFigureFormat;

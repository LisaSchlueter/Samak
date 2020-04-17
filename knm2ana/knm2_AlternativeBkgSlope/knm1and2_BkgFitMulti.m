% 
% First Tritium Background Slope Constraint
% Th. Lasserre, April 2020
%

% Read KNM1
% Read KNM1 data
range     = 5;
RunAnaArg = {'RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
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
    };%
NonPoissonScaleFactorKNM1 = 1.064;
KNM1                      = MultiRunAnalysis(RunAnaArg{:});
KNM1.exclDataStart        = KNM1.GetexclDataStart(range);
    
% Read KNM2 data
range         = 5;
RunAnaArg = {...
    'RunList','KNM2_Prompt',...        % all KNM2 golden runs
    'DataType','Real',...
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'ROIFlag','Default',...
    };%
NonPoissonScaleFactorKNM2 = 1.112;
KNM2                     = MultiRunAnalysis(RunAnaArg{:});
KNM2.exclDataStart       = KNM2.GetexclDataStart(range);

%%

% KNM2 data
nbin2  = numel(KNM2.RunData.qU(KNM2.exclDataStart:end));
hv2    = KNM2.RunData.qU(KNM2.exclDataStart:end);
c2     = KNM2.RunData.TBDIS(KNM2.exclDataStart:end)./KNM2.RunData.qUfrac(KNM2.exclDataStart:end)/KNM2.RunData.TimeSec;
ce2    = NonPoissonScaleFactorKNM2*KNM2.RunData.TBDISE(KNM2.exclDataStart:end)./KNM2.RunData.qUfrac(KNM2.exclDataStart:end)/KNM2.RunData.TimeSec;
% KNM1 data
nbin1  = numel(KNM1.RunData.qU(KNM1.exclDataStart:end));
hv1    = KNM1.RunData.qU(KNM1.exclDataStart:end);
c1     = KNM1.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1.RunData.TimeSec;
ce1    = NonPoissonScaleFactorKNM1*KNM1.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1.RunData.TimeSec;

% KNM1+2 data
hv     = [hv1 ; hv2];
c      = [c1 ; c2];
ce     = [ce1 ; ce2];
% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
Data = [hv/1e3 c*1e3 ce*1e3];
%Data = Data(1:end-1,:);

knm1b = 1:7;
%knm2b = 8:numel(Data(:,1));
knm2b = [8 14];

% Linear Model
parnames = {'baseline1','slope1','baseline2','slope2'}; 
ParIni   = [mean(Data(knm1b,2)) 0 mean(Data(knm2b,2)) 0]; 

% Fit1
Args = {ParIni, Data, '-c', 'min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2GaussKNM12',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline1 = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope1    = %.3f ± %.3f - limit knm1+2 < %.1f mcps/keV \n',par(2),err(2),norminv(0.68269,abs(par(:,2)),err(:,2)));
fprintf('  baseline2 = %.3f ± %.3f\n',par(3),err(3));
fprintf('  slope2   = %.3f ± %.3f - limit knm1+2 < %.1f mcps/keV \n',par(4),err(4),norminv(0.68269,abs(par(:,4)),err(:,4)));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,numel(knm1b)+numel(knm2b)-numel(ParIni));
fprintf('================================================================\n');

% Fit2
ParIni   = [par(1) 0 par(3) 0]; 
Args = {ParIni, Data, '-c', 'fix 1 3 ; min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2GaussKNM12',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline1 = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope1    = %.3f ± %.3f - limit knm1+2 < %.1f mcps/keV \n',par(2),err(2),norminv(0.68269,abs(par(:,2)),err(:,2)));
fprintf('  baseline2 = %.3f ± %.3f\n',par(3),err(3));
fprintf('  slope2   = %.3f ± %.3f - limit knm1+2 < %.1f mcps/keV \n',par(4),err(4),norminv(0.68269,abs(par(:,4)),err(:,4)));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,numel(knm1b)+numel(knm2b)-numel(ParIni));
fprintf('================================================================\n');


% Figure 1
close all
figure(1)
set(gcf, 'Position',  [100, 100, 1400, 1200])

subplot(2,1,1)
hdata = errorbar(Data(knm1b,1),Data(knm1b,2),Data(knm1b,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(knm1b,1),Model(par(1:2),Data(knm1b,1)),'LineWidth',3);
hold off
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'KNM1 Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV',par(1),err(1),par(2),err(2)),'Location','EastOutside')
grid on

subplot(2,1,2)
hdata = errorbar(Data(knm2b,1),Data(knm2b,2),Data(knm2b,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(knm2b,1),Model(par(3:4),Data(knm2b,1)),'LineWidth',3);
hold off
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],'KNM2 Data',sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV \n Upper Limit %0.2f mcps/keV',...
    par(3),err(3),par(4),err(4),norminv(0.68269,abs(par(:,4)),err(:,4))),'Location','EastOutside')
grid on

publish_figurePDF(1,'./knm1and2_BkgFitMulti.pdf');


%% Figure 2
figure(2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
corplot(errmat);
PrettyFigureFormat;

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
    'FSDFlag','Sibille0p5eV',...               % final state distribution
    'ELossFlag','KatrinT2',...            % energy loss function
    'SysBudget',22,...
    'DopplerEffectFlag','FSD_Knm1',...
    'Twin_SameCDFlag','ON',...
    'Twin_SameIsotopFlag','ON',...
    };%
NonPoissonScaleFactor = 1.064;
KNM1                     = MultiRunAnalysis(RunAnaArg{:});
KNM1.exclDataStart       = KNM1.GetexclDataStart(range);

% Read KNM2 data
range         = 5;
RunAnaArg = {...
    'RunList','KNM2_Prompt',...        % all KNM2 golden runs
    'DataType','Real',...
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'ROIFlag','Default',...
    };%
NonPoissonScaleFactor = 1.112;
KNM2                     = MultiRunAnalysis(RunAnaArg{:});
KNM2.exclDataStart       = KNM2.GetexclDataStart(range);

%%

% KNM2 data
nbin2  = numel(KNM2.RunData.qU(KNM2.exclDataStart:end));
hv2    = KNM2.RunData.qU(KNM2.exclDataStart:end);
c2     = KNM2.RunData.TBDIS(KNM2.exclDataStart:end)./KNM2.RunData.qUfrac(KNM2.exclDataStart:end)/KNM2.RunData.TimeSec;
ce2    = KNM2.RunData.TBDISE(KNM2.exclDataStart:end)./KNM2.RunData.qUfrac(KNM2.exclDataStart:end)/KNM2.RunData.TimeSec;
% KNM1 data
nbin1  = numel(KNM1.RunData.qU(KNM1.exclDataStart:end));
hv1    = KNM1.RunData.qU(KNM1.exclDataStart:end);
c1     = KNM1.RunData.TBDIS(KNM1.exclDataStart:end)./KNM1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1.RunData.TimeSec;
c1     = c1./mean(c1).*mean(c2);
%ce1    = KNM1.RunData.TBDISE(KNM1.exclDataStart:end)./KNM1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1.RunData.TimeSec;
ce1    = sqrt(KNM1.RunData.TBDIS(KNM1.exclDataStart:end)./mean(KNM1.RunData.TBDIS(KNM1.exclDataStart:end)).*mean(KNM1.RunData.TBDIS(KNM2.exclDataStart:end)))./KNM1.RunData.qUfrac(KNM1.exclDataStart:end)/KNM1.RunData.TimeSec;

% KNM1+2 data
hv     = [hv1 ; hv2];
c      = [c1 ; c2];
ce     = [ce1 ; ce2];
% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
Data = [hv/1e3 c*1e3 ce*1e3];
Data = Data(1:end-1,:);
% Linear Model
parnames = {'baseline','slope'}; 
ParIni   = [mean(Data(:,2)) 0]; 

% Case 1) Fit Overall Range
% Fit
Args                          = {ParIni, Data, '-c', 'min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope    = %.3f ± %.3f\n',par(2),err(2));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,numel(nbin1)+numel(nbin2)-numel(ParIni));
fprintf('================================================================\n');

% Figure 1
close all
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 500])
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks',...
       'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(:,1),Model(par,Data(:,1)),'LineWidth',3);
%hfiteu  =  plot(Data(:,1),Model(par-[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
%hfited  =  plot(Data(:,1),Model(par+[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
hold off;
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],...
    sprintf('KNM1+2 Data: %.0f counts in %.0f sec',...
    sum(KNM1.RunData.TBDIS(KNM1.exclDataStart:end)+KNM2.RunData.TBDIS(KNM2.exclDataStart:end)),...
    sum(KNM1.RunData.qUfrac(KNM1.exclDataStart:end).*KNM1.RunData.TimeSec+KNM2.RunData.qUfrac(KNM2.exclDataStart:end).*KNM2.RunData.TimeSec)),...
    sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV \n relative slope = %0.2f \\pm %0.2f %%/keV \n Upper Limit %0.2f mcps/keV',...
    par(1),err(1),par(2),err(2),par(2)/par(1)*100,err(2)/par(1)*100,norminv(0.68269,abs(par(:,2)),err(:,2))),'Location','SouthEast');
grid on
publish_figurePDF(1,'BKG_KNM1and2.pdf');

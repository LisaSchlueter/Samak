% Compare Samak - Fitrium response functions
%
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_systematics/results/'];
MakeDir(savedir);

%% load Fitrium
savenameF = sprintf('%ksn1_ResponseFunction_Fitrium.dat',savedir);

%% load Samak
savenameS = sprintf('%ksn1_ResponseFunction.mat',savedir);
if exist(savenameS,'file')
    dS = importdata(savenameS);
else
    RunAnaArg = {'RunList','KNM1',...
        'fixPar','E0 Norm Bkg',...
        'DataType','Real',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'TwinBias_Q',18573.73,...
        'SysBudget',24,...
        'pullFlag',99,...
        'NonPoissonScaleFactor',1};
    R = MultiRunAnalysis(RunAnaArg{:});
    
    
    qU = 18545;
    Te = (qU-5:0.1:qU+100)';
    RF = R.ModelObj.ComputeRF(Te,qU,'IntMode','Conv');
 
end


%% plot overlay
GetFigure;
pS = plot(Te-qU,RF,'-','LineWidth',2,'Color',rgb('Orange'));
hold on;
pF = plot(Te-qU,RF_fitrium,'-.','LineWidth',2,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('{\\itE}_{kin} - {\\itqU} (eV)'));
ylabel('Transmission probability');
legend('Samak','Fitrium','EdgeColor',rgb('Silver'));
%% Plot abs difference
GetFigure;
pS = plot(Te-qU,RF-RF_fitrium,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('{\\itE}_{kin} - {\\itqU} (eV)'));
ylabel('Transmission probability diff.');
legend('Samak - Fitrium','EdgeColor',rgb('Silver'));




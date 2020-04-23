% Investigate impact of new energy-loss parametrization
% twins with constant endpoint.

%% settings
ELossFlag1 = 'KatrinT2A20';
ELossFlag2 = 'KatrinT2';
range = 40; % eV below E0
savedir = [getenv('SamakPath'),'knm2ana/knm2_ElossKatrinT2A20/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_ElossComparison_%s_%s.mat',ELossFlag1,ELossFlag2);

if exist(savename,'file')
    load(savename)
else
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'fixPar','mNu E0 Bkg Norm',... % freePar
        'DataType','Twin',...
        'TwinBias_Q',18573.56,...
        'FSDFlag','BlindingKNM2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','ON',...
        'ISCSFlag','Edep'};
    
    %% init model
    E1 = MultiRunAnalysis(RunAnaArg{:},'ELossFlag',ELossFlag1);
    E1.exclDataStart = E1.GetexclDataStart(range);
    %% use MC data with elossflag1 as reference
    E1.InitModelObj_Norm_BKG;
    TBDIS_E1 = E1.ModelObj.TBDIS;
    E1.RunData.TBDIS = TBDIS_E1;
    E1.Fit;
    FitResults_ref = E1.FitResult;
    
    %%
    E1.ModelObj.ELossFlag = ELossFlag2;
    E1.ModelObj.InitializeELossFunction;
    E1.ModelObj.InitializeRF;
    E1.Fit;
    FitResults_El2 = E1.FitResult;
    
    save(savename,'FitResults_ref','FitResults_El2','ELossFlag1','ELossFlag2','RunAnaArg');
end

%% display results
fprintf('neutrino mass sq. shift = %.2g eV^2 (%s - %s) \n',...
    FitResults_ref.par(1)-FitResults_El2.par(1),ELossFlag1,ELossFlag2);
fprintf('endpoint shift          = %.2g eV   (%s - %s) \n',...
    FitResults_ref.par(2)-FitResults_El2.par(2),ELossFlag1,ELossFlag2);


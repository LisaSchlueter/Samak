% Test scattering transmission (aka detailed transmission)
% KNM2 uniform fit to 1 (random) single real run
% Lisa, Arpil 2020

%% settings
RunNr         = 56274; % 56274; 56275;
freePar       = 'E0 Bkg Norm';
DataType      = 'Real';
RecomputeFlag = 'ON';
%% load if possible
savedir = [getenv('SamakPath'),'knm2ana/knm2_ScatTF/results/'];
savename = sprintf('%sknm2_ScatTF_Run%.0f_%s_FitPar%s.mat',savedir,RunNr,DataType,strrep(freePar,' ',''));

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    RunAnaArg = {'RunNr',RunNr,...
        'fixPar',freePar,...     % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON'};
    
    % build object of MultiRunAnalysis class
    SR = RunAnalysis(RunAnaArg{:});
    SR.exclDataStart = SR.GetexclDataStart(40); % define fit range
    %% Using regular TF
    SR.ModelObj.AngularTFFlag = 'OFF';
    SR.ModelObj.recomputeRF='ON';
    SR.ModelObj.InitializeRF;
    SR.Fit;
    Te = SR.ModelObj.Te;
    qU = SR.ModelObj.qU;
    RFreg = SR.ModelObj.RF;
    FitResultReg = SR.FitResult;
    %% NEW: scattering TF
    SR.ModelObj.AngularTFFlag = 'ON';
    SR.ModelObj.recomputeRF='ON';
    SR.ModelObj.InitializeRF;
    SR.Fit;
    RFscat = SR.ModelObj.RF;
    FitResultScat = SR.FitResult;
    %% save
    save(savename,'FitResultReg','FitResultScat','Te','qU','RFreg','RFscat','RunAnaArg');
end

%% result
fprintf('E0 = %.3f eV (regular    TF) \n',FitResultReg.par(2)+18573.7);
fprintf('E0 = %.3f eV (scattering TF) \n',FitResultScat.par(2)+18573.7);
fprintf('E0 diff = %.3f eV (regular - scattering) \n',FitResultReg.par(2)-FitResultScat.par(2));


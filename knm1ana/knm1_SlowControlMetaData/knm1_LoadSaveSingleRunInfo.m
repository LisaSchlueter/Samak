% knm1 original settings
% load single run fit results (uniform, fixed nu-mass)
% load single run meta-data: column density, tritium purity, etc
% save for later use

savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file') && 1==2
    load(savefile);
else
    range = 40; % eV below the endpoint
    %% create MultiRunAnalysis object
    RunAnaArg =  {'RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','E0 Norm Bkg',...        % free Parameter!!
        'NonPoissonScaleFactor',1,...        % background uncertainty are (not) enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...          % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'DopplerEffectFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'DopplerEffectFlag','OFF'};   % already included inside FSD for knm1!!!!
    %%
    R = MultiRunAnalysis(RunAnaArg{:});
    R.exclDataStart = R.GetexclDataStart(range); % set region of interest

    FitResults = R.FitRunList;
    %%
    MakeDir(savedir);
    save(savefile,'FitResults','RunAnaArg','R');
end
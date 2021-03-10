RunLists = {'KNM1_part1_m149mvRW','KNM1_m149mvRW','KNM1_175mvRW','KNM1_m183mvRW','KNM1_300mvRW',...
    'KNM1_part1','KNM1_part2_1','KNM1_part2',...
    'KNM1-RhoD','KNM1-AntiRhoD',...
    'KNM1-TpurityLow','KNM1-TpurityHigh',...
    'KNM1upScan','KNM1downScan','KNM1rm',...
    'KNM1-FirstHalfTime','KNM1-MiddleHalfTime','KNM1-LastHalfTime',...
    'KNM1'};

fitresults = zeros(11,numel(RunLists));

for i=16:19
    M = MultiRunAnalysis('RunList',RunLists{i},...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','mNu E0 Norm Bkg eta',...        % free Parameter!!
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF');
    M.exclDataStart = M.GetexclDataStart(40);
    M.Fit;
    fitresults(1,i)  = M.FitResult.par(1);
    fitresults(2,i)  = M.FitResult.err(1);
    fitresults(3,i)  = M.ModelObj.Q_i+M.FitResult.par(2);
    fitresults(4,i)  = M.FitResult.err(2);
    fitresults(5,i)  = M.ModelObj.BKG_RateSec_i+M.FitResult.par(3);
    fitresults(6,i)  = M.FitResult.err(3);
    fitresults(7,i)  = M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1;
    fitresults(8,i)  = M.FitResult.err(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1);
    fitresults(9,i)  = M.FitResult.par(17).*1e10;
    fitresults(10,i) = M.FitResult.err(17).*1e10;
    fitresults(11,i) = M.FitResult.chi2min;
end

M = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','mNu E0 Norm Bkg eta',...        % free Parameter!!
        'AnaFlag','Ring',...
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF');
M.FitAllRings;
Nfit=2;
fitresults = zeros(11,Nfit);
for j=1:Nfit
    M = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Real',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','mNu E0 Norm Bkg eta',...    % free Parameter!!
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'fitter','minuit',...
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'Twin_SameCDFlag','OFF',...
        'Twin_SameIsotopFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF');

    M.exclDataStart = M.GetexclDataStart(40);
    M.ModelObj.Q_i=18573.73-10+20*rand;
    M.ModelObj.mnuSq_i=-2+4*rand;
    M.ModelObj.normFit_i=-0.02+0.04*rand;
    M.ModelObj.BKG_RateSec_i=0.295-0.02+0.04*rand;
    M.Fit;
    fitresults(1,j)  = M.FitResult.par(1);
    fitresults(2,j)  = M.FitResult.err(1);
    fitresults(3,j)  = M.ModelObj.Q_i+M.FitResult.par(2);
    fitresults(4,j)  = M.FitResult.err(2);
    fitresults(5,j)  = M.ModelObj.BKG_RateSec_i+M.FitResult.par(3);
    fitresults(6,j)  = M.FitResult.err(3);
    fitresults(7,j)  = M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1;
    fitresults(8,j)  = M.FitResult.err(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1);
    fitresults(9,j)  = M.FitResult.par(17).*1e10;
    fitresults(10,j) = M.FitResult.err(17).*1e10;
    fitresults(11,j) = M.FitResult.chi2min;
end
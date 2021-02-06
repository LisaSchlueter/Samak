initfile=@ref_RelicNuBkg_KNM1;
fitresults = zeros(11,1000);
for j=1
    M = RunAnalysis('RunNr',10,...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                'FakeInitFile',initfile,...
                'chi2','chi2Stat',...              % uncertainties: statistical or stat + systematic uncertainties
                'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                'fixPar','mNu E0 Norm Bkg eta',...    % free Parameter!!
                'NonPoissonScaleFactor',1,...     % background uncertainty are enhanced
                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                'FSDFlag','SibilleFull',...           % final state distribution
                'ELossFlag','KatrinT2',...            % energy loss function
                'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
                'DopplerEffectFlag','FSD',...
                'SynchrotronFlag','ON',...
                'AngularTFFlag','OFF');
        

    statfluct = zeros(numel(D.RunData.qU),1);
    for i=1:numel(D.RunData.qU)
        gm=gmdistribution(D.RunData.TBDIS(i),D.RunData.TBDIS(i));
        statfluct(i) = random(gm)-D.RunData.TBDIS(i);
    end
    M.RunData.TBDIS = M.RunData.TBDIS+statfluct;
    M.exclDataStart = M.GetexclDataStart(40);
    M.Fit;
    fitresults(j) = M.FitResult.par(17).*1e10;
end
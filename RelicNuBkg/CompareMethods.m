eta_Chi2 = zeros(1,20);
eta_fit  = zeros(1,20);
for i=1:21
    mNu = (i-1)/10;
    Q = RelicNuDebug('Params','KNM1');
    Q.Chi2Twin('Recompute','ON','Plot','OFF','Syst','ON','fitPar','mNu E0 Norm Bkg','DeltaChi2',1,'TwinBias_mnuSq',mNu);
    eta_Chi2(i) = Q.etaSensitivity;
    R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chiCMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','mNu E0 Norm Bkg eta',...        % free Parameter!!
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1,...     % background uncertainty are enhanced
        'fitter','minuit',...                 % minuit standard, matlab to be tried
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'Twin_SameCDFlag','OFF',...
        'Twin_SameIsotopFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'TwinBias_Q',18573.73,...
        'TwinBias_mnuSq',mNu);

    R.exclDataStart=R.GetexclDataStart(40);
    R.ModelObj.mnuSq_i = R.TwinBias_mnuSq;
    R.InitModelObj_Norm_BKG('Recompute','ON');
    R.Fit;
    R.FitResult
    eta_fit(i) = R.FitResult.err(17).*1e10;
end
save('./Results_mNuFree.mat','eta_fit','eta_Chi2');
mnu=linspace(0,2,21);
C=plot(mnu,eta_Chi2,'LineWidth',2);
hold on;
F=plot(mnu,eta_fit,'LineWidth',2);
xlabel('m_{\nu}^{2}');
ylabel('\eta');
legend('\chi^{2} Scan','\eta Fit');
legend boxoff;
PrettyFigureFormat;
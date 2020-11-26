M = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','mNu E0 Norm Bkg',...                   % free Parameter!!
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...          % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'Twin_SameCDFlag','OFF',...
        'Twin_SameIsotopFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'TwinBias_Q',18573.73,...
        'TwinBias_mnuSq',1);

M.exclDataStart = M.GetexclDataStart(40);
M.InitModelObj_Norm_BKG('Recompute','ON');

M.Fit;
ErrTotal = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('FSD','ON'),'BkgCM','OFF');
M.Fit;
ErrFSD = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('RF','ON'),'BkgCM','OFF');
M.Fit;
ErrRF = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('TASR','ON'),'BkgCM','OFF');
M.Fit;
ErrTASR = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
M.Fit;
ErrStack = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('FPDeff','ON'),'BkgCM','OFF');
M.Fit;
ErrFPD = M.FitResult.err(1);
M.ComputeCM('SysEffects',struct('TC','ON'),'BkgCM','OFF');
M.Fit;
ErrTC = M.FitResult.err(1);
M.ComputeCM('BkgCM','ON');
M.Fit;
ErrBkg = M.FitResult.err(1);
M.NonPoissonScaleFactor = 1;
M.ComputeCM('SysEffects',struct(),'BkgCM','OFF')
M.Fit;
ErrStat = M.FitResult.err(1);
ErrFSD = sqrt(ErrFSD.^2-ErrStat.^2);
ErrRF = sqrt(ErrRF.^2-ErrStat.^2);
ErrTASR = sqrt(ErrTASR.^2-ErrStat.^2);
ErrStack = sqrt(ErrStack.^2-ErrStat.^2);
ErrFPD = sqrt(ErrFPD.^2-ErrStat.^2);
ErrTC = sqrt(ErrTC.^2-ErrStat.^2);
ErrBkg = sqrt(ErrBkg.^2-ErrStat.^2);
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

ErrTotal = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('FSD','ON'),'BkgCM','OFF');
ErrFSD = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('RF','ON'),'BkgCM','OFF');
ErrRF = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('TASR','ON'),'BkgCM','OFF');
ErrTASR = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
ErrStack = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('FPDeff','ON'),'BkgCM','OFF');
ErrFPD = sqrt(diag(M.FitCM));
M.ComputeCM('SysEffects',struct('TC','ON'),'BkgCM','OFF');
ErrTC = sqrt(diag(M.FitCM));
M.ComputeCM('BkgCM','ON');
ErrBkg = sqrt(diag(M.FitCM));

T = plot(M.RunData.qU-M.ModelObj.Q,ErrTotal,'LineWidth',2);
hold on;
f = plot(M.RunData.qU-M.ModelObj.Q,ErrFSD,'LineWidth',2);
r = plot(M.RunData.qU-M.ModelObj.Q,ErrRF,'LineWidth',2);
t = plot(M.RunData.qU-M.ModelObj.Q,ErrTASR,'LineWidth',2);
s = plot(M.RunData.qU-M.ModelObj.Q,ErrStack,'LineWidth',2);
F = plot(M.RunData.qU-M.ModelObj.Q,ErrFPD,'LineWidth',2);
c = plot(M.RunData.qU-M.ModelObj.Q,ErrTC,'LineWidth',2);
b = plot(M.RunData.qU-M.ModelObj.Q,ErrBkg,'LineWidth',2);
l = legend([T f r t s F c b],'Total','FSD','RF','TASR','Stack','FPD','TC','Bkg');
PrettyFigureFormat;
hold off;
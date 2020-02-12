%% settings
RunList               = 'KNM1';
exclDataStart         = 14;
nTrials               = 5000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','OFF','FSD','OFF','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','OFF','TCoff_RAD','OFF','TCoff_OTHER','OFF');

%% Init Model Object and covariance matrix object
TwinFRON = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON');
TwinFRON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');

%% Display Effect
TwinFRON.FitCM_Obj.DisplayCMInfo

%% Plot Total Covariance Matrix
TwinFRON.FitCM_Obj.PlotCM('ConvergenceTest','OFF','qUWindowIndexMin',TwinFRON.exclDataStart,'qUWindowIndexMax',40,'Mode','Frac','saveplot','ON');

%% Plot Correlation Matrix
TwinFRON.PlotFitCovCorMatrices('Mode','CM');

%% Sensitivity - Via ASimov Fit Error
% Stat
TwinFRON.chi2='chi2Stat';TwinFRON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinFRON.Fit; TwinFRON.PlotFit;
Rstat = TwinFRON.FitResult;
% CM
TwinFRON.chi2='chi2CM';TwinFRON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinFRON.PlotFitCovCorMatrices('Mode','Frac');
TwinFRON.Fit; TwinFRON.PlotFit;
Rsys  = TwinFRON.FitResult;
% CM
TwinFRON.chi2='chi2CMShape';TwinFRON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinFRON.PlotFitCovCorMatrices('Mode','Shape');
TwinFRON.Fit; TwinFRON.PlotFit;
Rsysshape = TwinFRON.FitResult;

%%
t = PrintTable('TASR Sensitivity Fit Results');
t.addRow('Parameter','Value','Error','Value','Error','Value','Error');
t.addRow('m^2_\beta \, (eV^2)',...
    sprintf('%.2f',Rstat.par(1)),sprintf('%.2f',Rstat.err(1)),...
    sprintf('%.2f',Rsys.par(1)),sprintf('%.2f',Rsys.err(1)),...
    sprintf('%.2f',Rsysshape.par(1)),sprintf('%.2f',Rsysshape.err(1)));
t.addRow('E_{0,eff}  \, (eV)',...
    sprintf('%.2f',TwinFRON.ModelObj.Q_i+Rstat.par(2)),sprintf('%.2f',Rstat.err(2)),...
    sprintf('%.2f',TwinFRON.ModelObj.Q_i+Rsys.par(2)),sprintf('%.2f',Rsys.err(2)),...
    sprintf('%.2f',TwinFRON.ModelObj.Q_i+Rsysshape.par(2)),sprintf('%.2f',Rsysshape.err(2)));
t.addRow('Background  \, (mcps)',...
    sprintf('%.2f',(TwinFRON.ModelObj.BKG_RateSec_i + Rstat.par(3))*1e3),sprintf('%.2f',Rstat.err(3)*1e3),...
    sprintf('%.2f',(TwinFRON.ModelObj.BKG_RateSec_i + Rsys.par(3))*1e3),sprintf('%.2f',Rsys.err(3)*1e3),...
    sprintf('%.2f',(TwinFRON.ModelObj.BKG_RateSec_i + Rsysshape.par(3))*1e3),sprintf('%.2f',Rsysshape.err(3)*1e3));
t.addRow('Normalization',...
    sprintf('%.2f',Rstat.par(4)+1),sprintf('%.2f',Rstat.err(3)),...
    sprintf('%.2f',Rsys.par(4)+1),sprintf('%.2f',Rsys.err(3)),...
    sprintf('%.2f',Rsysshape.par(4)+1),sprintf('%.2f',Rsysshape.err(3)));
t.addRow('','','','','','','');
t.addRow('\chi^2/dof',...
    [sprintf('%.2f',Rstat.chi2min) '/' sprintf('%.2f',Rstat.dof)],'',...
    [sprintf('%.2f',Rsys.chi2min) '/' sprintf('%.2f',Rsys.dof)],'',...
    [sprintf('%.2f',Rsysshape.chi2min) '/' sprintf('%.2f',Rsysshape.dof)],'');
t.display;
t.HasHeader = true;
t.Format = 'tex';
t.Caption = sprintf('Fit Results  - Fitter = %s - Range = [%.1f - %.1f] eV',...
    TwinFRON.fitter,TwinFRON.ModelObj.qU(TwinFRON.exclDataStart),TwinFRON.ModelObj.qU(end));
t.print;
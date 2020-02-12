%% settings
RunList               = 'KNM1';
exclDataStart         = 2;
nTrials               = 1000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('FSD','OFF','RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','OFF','TCoff_RAD','ON','TCoff_OTHER','ON');

%% Init Model Object and covariance matrix object
TwinTCON = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON');
TwinTCON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');

%% Display Effect
TwinTCON.FitCM_Obj.DisplayCMInfo

%% Plot Total Covariance Matrix
TwinTCON.FitCM_Obj.PlotCM('ConvergenceTest','OFF','qUWindowIndexMin',TwinTCON.exclDataStart,'qUWindowIndexMax',40,'Mode','Frac','saveplot','ON');

%% Plot Correlation Matrix
TwinTCON.PlotFitCovCorMatrices('Mode','CM');

%% Sensitivity - Via ASimov Fit Error
% Stat
TwinTCON.chi2='chi2Stat';TwinTCON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinTCON.Fit; TwinTCON.PlotFit;
Rstat = TwinTCON.FitResult;
% CM
TwinTCON.chi2='chi2CM';TwinTCON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinTCON.PlotFitCovCorMatrices('Mode','Frac');
TwinTCON.Fit; TwinTCON.PlotFit;
Rsys  = TwinTCON.FitResult;
% CM
TwinTCON.chi2='chi2CMShape';TwinTCON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinTCON.PlotFitCovCorMatrices('Mode','Shape');
TwinTCON.Fit; TwinTCON.PlotFit;
Rsysshape = TwinTCON.FitResult;

%%
t = PrintTable('TASR Sensitivity Fit Results');
t.addRow('Parameter','Value','Error','Value','Error','Value','Error');
t.addRow('m^2_\beta \, (eV^2)',...
    sprintf('%.2f',Rstat.par(1)),sprintf('%.2f',Rstat.err(1)),...
    sprintf('%.2f',Rsys.par(1)),sprintf('%.2f',Rsys.err(1)),...
    sprintf('%.2f',Rsysshape.par(1)),sprintf('%.2f',Rsysshape.err(1)));
t.addRow('E_{0,eff}  \, (eV)',...
    sprintf('%.2f',TwinTCON.ModelObj.Q_i+Rstat.par(2)),sprintf('%.2f',Rstat.err(2)),...
    sprintf('%.2f',TwinTCON.ModelObj.Q_i+Rsys.par(2)),sprintf('%.2f',Rsys.err(2)),...
    sprintf('%.2f',TwinTCON.ModelObj.Q_i+Rsysshape.par(2)),sprintf('%.2f',Rsysshape.err(2)));
t.addRow('Background  \, (mcps)',...
    sprintf('%.2f',(TwinTCON.ModelObj.BKG_RateSec_i + Rstat.par(3))*1e3),sprintf('%.2f',Rstat.err(3)*1e3),...
    sprintf('%.2f',(TwinTCON.ModelObj.BKG_RateSec_i + Rsys.par(3))*1e3),sprintf('%.2f',Rsys.err(3)*1e3),...
    sprintf('%.2f',(TwinTCON.ModelObj.BKG_RateSec_i + Rsysshape.par(3))*1e3),sprintf('%.2f',Rsysshape.err(3)*1e3));
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
    TwinTCON.fitter,TwinTCON.ModelObj.qU(TwinTCON.exclDataStart),TwinTCON.ModelObj.qU(end));
t.print;
%% settings
RunList               = 'KNM1';
exclDataStart         = 2;
nTrials               = 5000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','ON','FSD','OFF','RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','OFF','TCoff_RAD','OFF','TCoff_OTHER','OFF');

%% Init Model Object and covariance matrix object
TwinRadON = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON');
TwinRadON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');

%% Display Effect
TwinRadON.FitCM_Obj.DisplayCMInfo

%% Plot Total Covariance Matrix
TwinRadON.FitCM_Obj.PlotCM('ConvergenceTest','OFF','qUWindowIndexMin',TwinRadON.exclDataStart,'qUWindowIndexMax',40,'Mode','Frac','saveplot','ON');

%% Plot Correlation Matrix
TwinRadON.PlotFitCovCorMatrices('Mode','CM');

%% Sensitivity - Via ASimov Fit Error
% Stat
TwinRadON.chi2='chi2Stat';TwinRadON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinRadON.Fit; TwinRadON.PlotFit;
Rstat = TwinRadON.FitResult;
% CM
TwinRadON.chi2='chi2CM';TwinRadON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinRadON.PlotFitCovCorMatrices('Mode','Frac');
TwinRadON.Fit; TwinRadON.PlotFit;
Rsys  = TwinRadON.FitResult;
% CM
TwinRadON.chi2='chi2CMShape';TwinRadON.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
TwinRadON.PlotFitCovCorMatrices('Mode','Shape');
TwinRadON.Fit; TwinRadON.PlotFit;
Rsysshape = TwinRadON.FitResult;

%%
t = PrintTable('TASR Sensitivity Fit Results');
t.addRow('Parameter','Value','Error','Value','Error','Value','Error');
t.addRow('m^2_\beta \, (eV^2)',...
    sprintf('%.2f',Rstat.par(1)),sprintf('%.2f',Rstat.err(1)),...
    sprintf('%.2f',Rsys.par(1)),sprintf('%.2f',Rsys.err(1)),...
    sprintf('%.2f',Rsysshape.par(1)),sprintf('%.2f',Rsysshape.err(1)));
t.addRow('E_{0,eff}  \, (eV)',...
    sprintf('%.2f',TwinRadON.ModelObj.Q_i+Rstat.par(2)),sprintf('%.2f',Rstat.err(2)),...
    sprintf('%.2f',TwinRadON.ModelObj.Q_i+Rsys.par(2)),sprintf('%.2f',Rsys.err(2)),...
    sprintf('%.2f',TwinRadON.ModelObj.Q_i+Rsysshape.par(2)),sprintf('%.2f',Rsysshape.err(2)));
t.addRow('Background  \, (mcps)',...
    sprintf('%.2f',(TwinRadON.ModelObj.BKG_RateSec_i + Rstat.par(3))*1e3),sprintf('%.2f',Rstat.err(3)*1e3),...
    sprintf('%.2f',(TwinRadON.ModelObj.BKG_RateSec_i + Rsys.par(3))*1e3),sprintf('%.2f',Rsys.err(3)*1e3),...
    sprintf('%.2f',(TwinRadON.ModelObj.BKG_RateSec_i + Rsysshape.par(3))*1e3),sprintf('%.2f',Rsysshape.err(3)*1e3));
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
    TwinRadON.fitter,TwinRadON.ModelObj.qU(TwinRadON.exclDataStart),TwinRadON.ModelObj.qU(end));
t.print;
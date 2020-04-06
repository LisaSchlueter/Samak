%
% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% 
% Last Updated: 15/08/2019
%

%% settings
RunList               = 'KNM1';
exclDataStart         = 13; % 27 subruns
RecomputeFlag         = 'OFF';
BkgCM                 = 'ON';

%% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar','mNu Norm E0 Bkg',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...%'SAENZ',...%Sibille0p5eV',...%'Sibille',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22);
% Modification - 24/10/2019
%    'fixPar','5 6 7 8 9 10 11',...
%% CM Shape
Real.chi2='chi2CMShape'; Real.ComputeCM('BkgCM',BkgCM);

%Real.PlotFitCovCorMatrices('Mode','Shape');
Real.Fit('CATS','OFF');
Real.PlotFit('LabelFlag','data','saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
%Real.PlotFit('LabelFlag','FinalKNM1','saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
%Real.PlotFit('LabelFlag','FinalKNM1','saveplot','pdf','ErrorBarScaling',50,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
return;

%Real.PlotFit('LabelFlag','FinalKNM1','saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','BW');
%Real.PlotFit('LabelFlag','FinalKNM1','saveplot','pdf','ErrorBarScaling',50,'YLimRes',[-2.2,2.9],'Colors','BW');
Rsysshape = Real.FitResult;
Real.PlotFitCovCorMatrices('Mode','Shape');

%% Sensitivity - Via ASimov Fit Error
% Stat
Real.chi2='chi2Stat'; Real.ComputeCM('BkgCM',BkgCM);
Real.Fit;
Rstat = Real.FitResult;

%% CM
SysEffects = Real.GetDefaultEffects();
Real.chi2='chi2CM';Real.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
Real.PlotFitCovCorMatrices('Mode','Frac');
Real.Fit; Real.PlotFit('saveplot','png');
Rsys  = Real.FitResult;


%% Plot Correlation EO/m2
Real.PlotFitm2e0Contours('ErrorOnModel','OFF');

%%
t = PrintTable('Sensitivity Fit Results');
t.addRow('Parameter','Value','Error','Value','Error','Value','Error');
t.addRow('m^2_\beta \, (eV^2)',...
    sprintf('%.2f',Rstat.par(1)),sprintf('%.2f',Rstat.err(1)),...
    sprintf('%.2f',Rsys.par(1)),sprintf('%.2f',Rsys.err(1)),...
    sprintf('%.2f',Rsysshape.par(1)),sprintf('%.2f',Rsysshape.err(1)));
t.addRow('E_{0,eff}  \, (eV)',...
    sprintf('%.2f',Real.ModelObj.Q_i+Rstat.par(2)),sprintf('%.2f',Rstat.err(2)),...
    sprintf('%.2f',Real.ModelObj.Q_i+Rsys.par(2)),sprintf('%.2f',Rsys.err(2)),...
    sprintf('%.2f',Real.ModelObj.Q_i+Rsysshape.par(2)),sprintf('%.2f',Rsysshape.err(2)));
t.addRow('\rm{Background}  \, (mcps)',...
    sprintf('%.2f',(Real.ModelObj.BKG_RateSec_i + Rstat.par(3))*1e3),sprintf('%.2f',Rstat.err(3)*1e3),...
    sprintf('%.2f',(Real.ModelObj.BKG_RateSec_i + Rsys.par(3))*1e3),sprintf('%.2f',Rsys.err(3)*1e3),...
    sprintf('%.2f',(Real.ModelObj.BKG_RateSec_i + Rsysshape.par(3))*1e3),sprintf('%.2f',Rsysshape.err(3)*1e3));
t.addRow('\rm{Normalization}',...
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
    Real.fitter,Real.ModelObj.qU(Real.exclDataStart),Real.ModelObj.qU(end));
t.print;
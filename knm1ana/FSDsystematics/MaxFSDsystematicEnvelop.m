%% settings
RunList               = 'KNM1';
exclDataStart         = 2;
nTrials               = 5000;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','OFF','FSD','ON','RF_RX','OFF','RF_EL','OFF','RF_BF','OFF','BkgShape','OFF','TCoff_RAD','OFF','TCoff_OTHER','OFF');

%% Init Model Object and covariance matrix object
if ~exist('TwinFSD','var')
TwinFSD = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON');
end

%% FSD Blinding
countfit=0;
TwinFSD.SimulateStackRuns('FSDflag',0);
for exclDataStart_i = exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU)-15)
    countfit=countfit+1;
    TwinFSD.exclDataStart = exclDataStart_i;
      switch  TwinFSD.chi2
        case {'chi2CM','chi2CMShape'}
            TwinFSD.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
    end
    TwinFSD.Fit;
    FSDResults(countfit) = TwinFSD.FitResult;
end
% Format Results
mass2_fsdB    = cell2mat(arrayfun(@(x) x.par(1),FSDResults,'UniformOutput',false));
e0_fsdB       = cell2mat(arrayfun(@(x) x.par(2) + TwinFSD.ModelObj.Q_i,FSDResults,'UniformOutput',false));
B_fsdB        = cell2mat(arrayfun(@(x) x.par(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
N_fsdB        = cell2mat(arrayfun(@(x) x.par(4),FSDResults,'UniformOutput',false));
mass2Err_fsdB = cell2mat(arrayfun(@(x) x.err(1),FSDResults,'UniformOutput',false));
e0Err_fsdB    = cell2mat(arrayfun(@(x) x.err(2) ,FSDResults,'UniformOutput',false));
BErr_fsdB     = cell2mat(arrayfun(@(x) x.err(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
NErr_fsdB     = cell2mat(arrayfun(@(x) x.err(4),FSDResults,'UniformOutput',false));
chi2_fsdB     = cell2mat(arrayfun(@(x) x.chi2min,FSDResults,'UniformOutput',false));

%% FSD FT/SAENZ
countfit=0;
TwinFSD.SimulateStackRuns('FSDflag',1);
for exclDataStart_i = exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU)-15)
    countfit=countfit+1;
    TwinFSD.exclDataStart = exclDataStart_i;
    switch  TwinFSD.chi2
        case {'chi2CM','chi2CMShape'}
            TwinFSD.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
    end
    TwinFSD.Fit;
    FSDResults(countfit) = TwinFSD.FitResult;
end
% Format Results
mass2_fsdFT    = cell2mat(arrayfun(@(x) x.par(1),FSDResults,'UniformOutput',false));
e0_fsdFT       = cell2mat(arrayfun(@(x) x.par(2) + TwinFSD.ModelObj.Q_i,FSDResults,'UniformOutput',false));
B_fsdFT        = cell2mat(arrayfun(@(x) x.par(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
N_fsdFT        = cell2mat(arrayfun(@(x) x.par(4),FSDResults,'UniformOutput',false));
mass2Err_fsdFT = cell2mat(arrayfun(@(x) x.err(1),FSDResults,'UniformOutput',false));
e0Err_fsdFT    = cell2mat(arrayfun(@(x) x.err(2) ,FSDResults,'UniformOutput',false));
BErr_fsdFT     = cell2mat(arrayfun(@(x) x.err(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
NErr_fsdFT     = cell2mat(arrayfun(@(x) x.err(4),FSDResults,'UniformOutput',false));
chi2_fsdFT     = cell2mat(arrayfun(@(x) x.chi2min,FSDResults,'UniformOutput',false));

%% FSD DOSS
countfit=0;
TwinFSD.SimulateStackRuns('FSDflag',2);
for exclDataStart_i = exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU)-15)
    countfit=countfit+1;
    TwinFSD.exclDataStart = exclDataStart_i;
      switch  TwinFSD.chi2
        case {'chi2CM','chi2CMShape'}
            TwinFSD.ComputeCM('SysEffects',SysEffects,'BkgCM','OFF');
    end
    TwinFSD.Fit;
    FSDResults(countfit) = TwinFSD.FitResult;
end
% Format Results
mass2_fsdD    = cell2mat(arrayfun(@(x) x.par(1),FSDResults,'UniformOutput',false));
e0_fsdD       = cell2mat(arrayfun(@(x) x.par(2) + TwinFSD.ModelObj.Q_i,FSDResults,'UniformOutput',false));
B_fsdD        = cell2mat(arrayfun(@(x) x.par(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
N_fsdD        = cell2mat(arrayfun(@(x) x.par(4),FSDResults,'UniformOutput',false));
mass2Err_fsdD = cell2mat(arrayfun(@(x) x.err(1),FSDResults,'UniformOutput',false));
e0Err_fsdD    = cell2mat(arrayfun(@(x) x.err(2) ,FSDResults,'UniformOutput',false));
BErr_fsdD     = cell2mat(arrayfun(@(x) x.err(3) + TwinFSD.ModelObj.BKG_RateAllFPDSec,FSDResults,'UniformOutput',false));
NErr_fsdD     = cell2mat(arrayfun(@(x) x.err(4),FSDResults,'UniformOutput',false));
chi2_fsdD     = cell2mat(arrayfun(@(x) x.chi2min,FSDResults,'UniformOutput',false));

%% Comparison - numasssquared shift verus range for the 3 cases
qUmin = round(TwinFSD.ModelObj.qU(TwinFSD.exclDataStart));
qUmax = round(TwinFSD.ModelObj.qU(end));
range    = round(TwinFSD.ModelObj.qU(TwinFSD.exclDataStart)-TwinFSD.ModelObj.Q_i);
myMainTitle=[sprintf('KATRIN - KNM1 Uniform Stacked Twin Toy Monte Carlo - %d Runs - [%.0f - %0.f] eV - Stat. Only',...
    numel(TwinFSD.RunList),qUmin,qUmax)];
maintitle=myMainTitle;
savefile=sprintf('plots/KNM1_TwinStackedToyMC_FSDbias_%d_%s_%.0feVbelowE0-1.png',...
    numel(TwinFSD.RunList),TwinFSD.chi2,abs(range));
fig1 = figure('Name','KATRIN - KNM1 Uniform Stacked Fit','NumberTitle','off','rend','painters','pos',[10 10 1400 700]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

l1=stairs(TwinFSD.RunData.qU(exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU))-15)-TwinFSD.ModelObj.Q_i,mass2_fsdB(1:12),'LineWidth',5);
hold on
l2=stairs(TwinFSD.RunData.qU(exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU))-15)-TwinFSD.ModelObj.Q_i,mass2_fsdFT(1:12)*randn(1),'LineWidth',5);
l3=stairs(TwinFSD.RunData.qU(exclDataStart:2:numel(exclDataStart:numel(TwinFSD.RunData.qU))-15)-TwinFSD.ModelObj.Q_i,mass2_fsdD(1:12)*randn(1),'LineWidth',5);
hold off
a = legend([l1 l2 l3],'Data=Model=Blinded FSD','FSD SAENZ','FSD DOSS','Location','NorthEast');
            legend('boxoff');
            xlabel(sprintf('retarding energy - %.1 eVf',TwinFSD.ModelObj.Q_i));
ylabel('\delta m_\beta^2 x random number (eV^2)')
PrettyFigureFormat; set(gca,'FontSize',24);
export_fig(savefile,'-png','-r300');

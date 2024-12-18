% KNM1 fit qU-background-slope with cut offs

%% settings
RecomputeFlag = 'OFF';
Plots        = 'OFF';
DataType     = 'Real';
range        = 40;
freePar      = 'mNu E0 Norm Bkg';
savedir = [getenv('SamakPath'),'knm1ana/knm1_FitBkgSlope/results/'];
MaxSlopeCpsPereV = 15*1e-06;  %cut-off for linear fit (cov-mat)
savefile = sprintf('%sknm1_%s_%s_%.0feV_%.1fmuCpsPerS.mat',savedir,DataType,strrep(freePar,' ',''),range,1e6.*MaxSlopeCpsPereV);

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    
    %% Init Model Object and covariance matrix object
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',1,...
        'SysBudget',22,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'AngularTFFlag','OFF',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','OFF');
    
    Real.exclDataStart =    Real.GetexclDataStart(range);
    Real.chi2 = 'chi2Stat';
    Real.InitModelObj_Norm_BKG;
    Real.Fit;
    FitResultStat = Real.FitResult;
    
    Real.chi2 = 'chi2CMShape';
    Real.ComputeCM('SysEffects',struct('FSD','OFF'),...
        'BkgPtCM','OFF',...
        'BkgCM','ON',...
        'BkgMode','SlopeFit',...
        'MaxSlopeCpsPereV',99,...%MaxSlopeCpsPereV,...
        'BkgRange',-5,...
        'PlotSaveCM','ON'); 
    Real.Fit;
    FitCM = Real.FitCMShape;
    FitCMFrac = Real.FitCMFracShape;
    FitResult = Real.FitResult;
    mNuSq_sys = sqrt(FitResult.err(1).^2-FitResultStat.err(1).^2);
    
    save(savefile,'mNuSq_sys','MaxSlopeCpsPereV','FitResultStat','FitResult','Real');
end



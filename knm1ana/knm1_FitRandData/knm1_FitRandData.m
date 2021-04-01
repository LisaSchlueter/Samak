%
% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% Last Updated: 23/03/2021

%% settings
nFit = 2;
RecomputeFlag = 'OFF';
Plots        = 'OFF';
DataType     = 'Twin';
range        = 40;
chi2         = 'chi2CMShape';
freePar      = 'mNu E0 Norm Bkg';
FSDFlag      = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;

SysBudget    = 24;
savedir = [getenv('SamakPath'),'knm1ana/knm1_FitRandData/results/'];
TwinmNuSq = -0.97;
savefile = sprintf('%sknm1FitRandData_%s_%s_NP%.4g_%s_%.0feV_%s_TwinmNuSq%.3geV2_%.0ffit.mat',...
    savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag,TwinmNuSq,nFit);

if strcmp(chi2,'chi2CMShape')
    savefile = strrep(savefile,chi2,sprintf('%s_SysBudget%.0f',chi2,SysBudget));
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    
    %% Init Model Object and covariance matrix object
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'SysBudget',SysBudget,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2',...
        'AngularTFFlag','OFF',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','OFF',...
        'TwinBias_mnuSq',TwinmNuSq);
    
    Real.exclDataStart =    Real.GetexclDataStart(range);
    
    if strcmp(chi2,'chi2Stat')
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    else
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        Real.ComputeCM;
    end
    
    [FitPar, FitErr, FitChi2min, dof,TBDIS]  = Real.FitTwin('nSamples',nFit);
    
    MakeDir(savedir);
    save(savefile,'FitResult','Real','FitPar', 'FitErr', 'FitChi2min', 'dof','TBDIS');
end

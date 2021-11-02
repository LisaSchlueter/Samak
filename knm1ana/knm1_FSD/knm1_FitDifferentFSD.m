% Uniform fit using different FSD types
% Uniform Fit
% Golden Run List
% Golden Pixel List
% Last Updated: 6/10/2021

%% settings
RecomputeFlag = 'OFF';
DataType     = 'Real';
range        = 40;
chi2         = 'chi2Stat';
freePar      = 'mNu E0 Norm Bkg';
FSDFlag      = 'Sibille0p5eV';% KNM2';%
NPFactor = 1.064;
savedir = [getenv('SamakPath'),'knm1ana/knm1_FSD/results/'];
savefile = sprintf('%sknm1_FitDifferentFSD_%s_%s_NP%.4g_%s_%.0feV_%s.mat',...
    savedir,DataType,chi2,NPFactor,strrep(freePar,' ',''),range,FSDFlag);

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    
    if contains(FSDFlag,'Sibille')
        DopplerEffectFlag = 'OFF';
    else
        DopplerEffectFlag = 'FSD';
    end
    
    %% Init Model Object and covariance matrix object
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NPFactor,...
        'SysBudget',SysBudget,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2',...
        'AngularTFFlag','OFF',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag',DopplerEffectFlag);
    
   
    Real.exclDataStart =    Real.GetexclDataStart(range);
    
    if strcmp(chi2,'chi2Stat')
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    else
        Real.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        Real.ComputeCM;
    end
    
    Real.Fit;
    
    FitResult = Real.FitResult;
    save(savefile,'FitResult','Real');
end
%%
fprintf('m^2 = %.3f (%.3f +%.3f) eV^2 \n',FitResult.par(1),FitResult.errNeg(1),FitResult.errPos(1))
fprintf('E0 = %.3f  (%.3f +%.3f) eV \n',FitResult.par(2)+Real.ModelObj.Q_i,FitResult.errNeg(2),FitResult.errPos(2))
fprintf('B  = %.1f (%.1f +%.1f) mcps \n',(FitResult.par(3)+Real.ModelObj.BKG_RateSec_i)*1e3,1e3*FitResult.errNeg(3),1e3*FitResult.errPos(3))
fprintf('chi^2 = %.1f (%.0f dof) \n',FitResult.chi2min,FitResult.dof)



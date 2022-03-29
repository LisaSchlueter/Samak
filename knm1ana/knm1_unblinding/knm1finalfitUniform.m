%
% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List
% Last Updated: 23/03/2021

%% settings
RecomputeFlag = 'OFF';
Plots        = 'OFF';
DataType     = 'Real';
range        = 40;
chi2         = 'chi2CMShape';
freePar      = 'mNu E0 Norm Bkg';
FSDFlag      = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
SysBudget    = 22;
savedir = [getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'];
savefile = sprintf('%sknm1finalfitUniform_%s_%s_NP%.4g_%s_%.0feV_%s.mat',savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag);

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
        'DopplerEffectFlag','OFF');
    
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
fprintf('m^2 = %.3f (+-%.3f %.3f +%.3f) eV^2 \n',FitResult.par(1),(FitResult.errPos(1)-FitResult.errNeg(1))/2,FitResult.errNeg(1),FitResult.errPos(1))
fprintf('E0 = %.3f  (+-%.3f) eV \n',FitResult.par(2)+Real.ModelObj.Q_i,FitResult.err(2))
fprintf('B  = %.1f (%.1f +%.1f) mcps \n',(FitResult.par(3)+Real.ModelObj.BKG_RateSec_i)*1e3,1e3*FitResult.errNeg(3),1e3*FitResult.errPos(3))
fprintf('chi^2 = %.1f (%.0f dof) \n',FitResult.chi2min,FitResult.dof)
if strcmp(Plots,'ON')
    Real.PlotFit('LabelFlag','data','saveplot','pdf','ErrorBarScaling',1,'YLimRes',[-2.2,2.9],'Colors','RGB','DisplayStyle','PRL');
end


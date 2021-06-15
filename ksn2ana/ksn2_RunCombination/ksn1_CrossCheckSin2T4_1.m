% calculate chi^2 for sin2t4 = 1 and m4Sq = 0.28 (KSN-2 best fit)
mnu4Sq =0.2761;
sin2T4 = 1;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir)
savefile = sprintf('%sksn1_CrossCheckSin2T4_%.3g_m4Sq%.1feV2.mat',savedir,sin2T4,mnu4Sq);

if exist(savefile,'file')
    load(savefile)
else
    %% settings that might change
    nGridSteps = 30;
    DataType              = 'Real';
    range                 = 40;
    chi2                  = 'chi2CMShape';
    freePar               = 'E0 Norm Bkg';
    
    % configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.064;
    end
    Real = MultiRunAnalysis('RunList','KNM1',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'SysBudget',200,...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','FSD',...
        'BKG_PtSlope',-2.2*1e-06);
    Real.exclDataStart = Real.GetexclDataStart(range);
    
    Real.ModelObj.SetFitBiasSterile(mnu4Sq,sin2T4)
    Real.Fit;
    
    FitResult = Real.FitResult;
    chi2 = FitResult.chi2min;
    save(savefile,'FitResult','chi2')
    fprintf('save result to %s \n',savefile);
end

fprintf('fixed nu-mass \n')
fprintf('KSN-1 chi2min  = %.2f (m4^2 = %.2f sin2t4 = %.3g)       -> ksn2 minimum\n',chi2,mnu4Sq,sin2T4)
fprintf('KSN-1 chi2min  = %.2f (m4^2 = %.2f sin2t4 = %.3g)  -> actual minimum \n',21.4,77.1,0.031)

fprintf('Common chi2min = %.2f (m4^2 = %.2f sin2t4 = %.3g)  -> actual minimum \n',50.4,59.5,0.011)
fprintf('Common chi2min = %.2f (m4^2 = %.2f sin2t4 = %.3g)       -> ksn-2 minimum \n',chi2+27.5,0.28,1)

% fixed nu-mass

chi2min_ksn12_real = 50.4; % actual minimum of common chi^2 map

chi2min_ksn2 = 27.5;       % compare to common chi2 at m4^2 = 0.3 eV^2 sin2t4 = 1.000 
chi2min_ksn12_k2bf = chi2min_ksn2+chi2;

% best fit ksn-1 re-analysis: chi2min = 21.4  m4^2=77.1 eV^2    sin2T4 = 0.031
% compare to ksn-2 best fit : chi2min = 27.5  m4^2 = 0.3 eV^2   sin2t4 = 1.000 

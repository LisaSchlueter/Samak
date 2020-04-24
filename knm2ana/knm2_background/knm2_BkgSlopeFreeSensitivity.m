%% sensitivity on background slope from KNM2 data
%% settings
range = 40;
NonPoissonScaleFactor = 1.112;

%% load or calculate
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_BkgSlopeFreeSensitivity_NPfactor%.2f_%.0feV.mat',savedir,NonPoissonScaleFactor,range);
if exist(savename,'file')
    load(savename)
else
    %% set up model
    RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','mNu E0 Bkg Norm',...
    'DataType','Twin',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','ON',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.56,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor};
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    %% fit with free nu-mass + fixed bkg slope to 0
    T.fixPar = 'mNu E0 Bkg Norm';
    T.InitFitPar;
    T.Fit;
    FitResults_mNuSqFree = T.FitResult;
    
    %% fit with free nu-mass + free bkg slope
    T.fixPar = 'mNu E0 Bkg Norm BkgSlope';
    T.InitFitPar;
    T.Fit;
    FitResults_mNuSqBkgFree = T.FitResult;
    
    %% fit with fixed nu-mass + free bkg slope
    T.fixPar = 'E0 Bkg Norm BkgSlope';
    T.InitFitPar;
    T.Fit;
    FitResults_BkgFree = T.FitResult;

    save(savename,'FitResults_BkgFree','FitResults_mNuSqBkgFree','range',...
                  'FitResults_mNuSqFree','RunAnaArg','range','NonPoissonScaleFactor')
end

%% print results
fprintf('-----------------------------------------------------------\n')
fprintf('Background slope -> free mNuSq  = %.2f +- %.1f mcps / keV \n',FitResults_mNuSqBkgFree.par(12)*1e6,FitResults_mNuSqBkgFree.err(12)*1e6);
fprintf('Background slope -> fixed mNuSq =  %.2f +- %.1f mcps / keV \n',FitResults_BkgFree.par(12)*1e6,FitResults_BkgFree.err(12)*1e6);
fprintf('-----------------------------------------------------------\n')
fprintf('Neutrino mass sq.->  free bkg slope = %.3f +- %.3f eV^2 \n',FitResults_mNuSqBkgFree.par(1),0.5*(-FitResults_mNuSqBkgFree.errNeg(1)+FitResults_mNuSqBkgFree.errPos(1)))
fprintf('Neutrino mass sq.-> fixed bkg slope = %.3f +- %.3f eV^2 \n',FitResults_mNuSqFree.par(1),0.5*(-FitResults_mNuSqFree.errNeg(1)+FitResults_mNuSqFree.errPos(1)))
fprintf('-----------------------------------------------------------\n')
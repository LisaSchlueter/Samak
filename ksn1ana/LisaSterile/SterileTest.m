% test parallel grid search for sterile analysis
% Lisa, April 2020
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
MakeDir(savedir);
range = 40;
nGridSteps = 50;
NonPoissonScaleFactor=1;
chi2 = 'chi2Stat';
DataType = 'Real';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
savefile = sprintf('%sSterileTest_%s_%s_%s_%.0feVrange_%.0fGridSteps.mat',...
    savedir,RunList,DataType,strrep(freePar,' ',''),range,nGridSteps);
if exist(savefile,'file')
else
    RunAnaArg = {'RunList',RunList,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor};
    tic;
    T = MultiRunAnalysis(RunAnaArg{:});
    toc;
    T.exclDataStart = T.GetexclDataStart(range);
    
    T.Fit;
    FitResults_ref = T.FitResult;
    chi2_ref       = T.FitResult.chi2min;
    %% define grid
    mnu4Sq = repmat(linspace(0,50^2,nGridSteps)',1,nGridSteps);
    sin2T4 = repmat(linspace(0.001,1,nGridSteps),nGridSteps,1); %logspace(-3,0,nGridSteps)
    
    % make copy of models for parallel computing
    D = copy(repmat(T,nGridSteps,nGridSteps));
    arrayfun(@(x) x.SimulateStackRuns,D'); % get model object
    arrayfun(@(x,mu,sin) x.ModelObj.SetFitBiasSterile(mu,sin),D',mnu4Sq,sin2T4); % set sterile parameters
    
    %%  start timed grid search
    tic;
    arrayfun(@(x) x.Fit,D)
    toc;
    %%
    chi2Grid        = arrayfun(@(x) x.FitResult.chi2min,D);
    FitResultsGrid  = arrayfun(@(x) x.FitResult,D);
    save(savefile,'chi2_ref','FitResults_ref','RunAnaArg',...
        'chi2Grid','mnu4Sq','sin2T4','FitResultsGrid');
    fprintf('save file to %s \n',savefile)
end



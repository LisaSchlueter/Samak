file=load('FitRunListResults_KNM1_m149mvRW_chi2Stat_all_94eVbelowE0_fixPar15678910.mat');

counter=0;
for run=1:1:numel(file.Runs)
    counter=counter+1;
    Run      = (file.Runs(counter));
    mnuSq    = (file.FitResults.mnuSq(counter));
    mnuSqErr = (file.FitResults.mnuSqErr(counter));
    E0       = (file.FitResults.E0(counter));
    E0Err    = (file.FitResults.E0Err(counter));
    B        = (file.FitResults.B(counter));
    BErr     = (file.FitResults.BErr(counter));
    N        = (file.FitResults.N(counter));
    NErr     = (file.FitResults.NErr(counter));
    chi2min  = (file.FitResults.chi2min(counter));
    dof      = (file.FitResults.dof(counter));
    pValue   = (file.FitResults.pValue(counter));

    folderO  = '../../../tritium-data/fit/';  
    filename = [folderO num2str(Run) '_chi2Stat_94bE0_fixPar15678910.mat'];
   
    save(filename,'Run','mnuSq','mnuSqErr','E0','E0Err',...
        'B','BErr','N','NErr','chi2min','dof','pValue');
    disp(['Run Wise Fit Resulted Saved in ' filename]);
   
end
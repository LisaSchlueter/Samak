% get nu-mass sensitivity with fits with free nu-mass

chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';
SavePlt = 'ON';

nGridSteps = 2;%10;
sin2T4_i = logspace(-3,log10(0.5),nGridSteps);
mNu4Sq_i = logspace(-1,log10(40^2),nGridSteps);

PullFlag = 27;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_NuMassSensitivityFits_freemNuSq_%s_Pull%.0f_%s_nGrid%.0f.mat',...
    savedir,DataType,PullFlag,chi2,nGridSteps);

if exist(savefile,'file')
    load(savefile,'mNuSq','chi2min');
else
 
    %% configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg sin2T4 mnu4Sq',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',PullFlag,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(40);
    
    FitResults = cell(nGridSteps,nGridSteps);
    chi2min    = zeros(nGridSteps,nGridSteps);
   
    progressbar('Nu-mass sensitivity...');
    
    
    for i = 1:nGridSteps
        progressbar((i-1)./nGridSteps);
        for j=1:nGridSteps
            A.ModelObj.SetFitBias(0); % re-set fit bias
            A.ModelObj.SetFitBiasSterile(mNu4Sq_i(j),sin2T4_i(i)); % init sterile parameters
            A.Fit;
            FitResults{i,j} = A.FitResult;
            chi2min(i,j) = A.FitResult.chi2min;  
        end
    end
    
    save(savefile,'FitResults','chi2min','RunAnaArg','sin2T4_i','mNu4Sq_i','nGridSteps'); 
end

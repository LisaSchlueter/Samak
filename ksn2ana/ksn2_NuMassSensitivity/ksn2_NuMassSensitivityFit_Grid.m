% get nu-mass sensitivity with fits -> chi^2 profile a la grid search


chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';
SavePlt = 'ON';
PullFlag = 27;
nGridSteps = 10;

if strcmp(DataType,'Real')
    mNuSq = -1:0.2:2;%.5;
else
    mNuSq = -1:0.1:1;
end


mNu4Sq_i = logspace(-1,log10(40^2),nGridSteps);
sin2T4_i = logspace(-3,log10(0.5),nGridSteps);

%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end

RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg sin2T4 mnu4Sq',...%free par
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

mNu4Sq_i      = repmat(mNu4Sq_i,nGridSteps,1);
sin2T4_i      = repmat(sin2T4_i',1,nGridSteps);

mNu4Sq_grid    = reshape(mNu4Sq_i,nGridSteps^2,1);
sin2T4_grid    = reshape(sin2T4_i,nGridSteps^2,1);

%% make copy of model for parallel computing
D = copy(repmat(A,nGridSteps^2,1));
D = reshape(D,numel(D),1);


savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_NuMassSensitivityFits_%s_Pull%.0f_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f_Grid%.0f.mat',...
    savedir,DataType,PullFlag,chi2,min(mNuSq),max(mNuSq),numel(mNuSq),nGridSteps);

FitResults = cell(numel(mNuSq),nGridSteps^2);
chi2min    = zeros(numel(mNuSq),nGridSteps^2);


parfor j =1:nGridSteps^2
    mNu4Sq_par = mNu4Sq_grid(j); % 95;%logspace(-1,log10(40^2),nGridSteps);
    sin2T4_par = sin2T4_grid(j);%0.02;%logspace(-3,log10(0.5),nGridSteps);

    progressbar('Nu-mass sensitivity...');
    for i = 1:numel(mNuSq)
        progressbar((i-1)./numel(mNuSq));
        D(j).ModelObj.SetFitBias(0); % reset
        D(j).ModelObj.SetFitBiasSterile(mNu4Sq_par,sin2T4_par);
        D(j).ModelObj.mnuSq_i = mNuSq(i);
        D(j).Fit;
        FitResults{i,j} = D(j).FitResult;
        chi2min(i,j) = D(j).FitResult.chi2min;
    end
   
end
 save(savefile,'FitResults','mNuSq','chi2min','RunAnaArg','mNu4Sq_i','sin2T4_i');

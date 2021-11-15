

% compute covariance matrix for varyying FSD, based on trueness

% get correlation matrix
 nSamples = 1e3;
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefileResult = sprintf('%sknm1_GasComposition_Sampling%.0f.mat',savedir,nSamples);
if exist(savefileResult,'file') 
    dR = importdata(savefileResult);
else
    
    % Init Model Object and covariance matrix object
    R = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF');
    R.exclDataStart = R.GetexclDataStart(40);  % 40eV range = 27 subruns
 
    %% get covariance matrix
    savefileCovMat = sprintf('%sknm1_LARADataCovMat_5000samples.mat',savedir);
    dCov = importdata(savefileCovMat);
    
   
    mNuSq = zeros(nSamples,1);
    FitResults = cell(nSamples,1);
    for i=1:nSamples
        progressbar(i/nSamples)
        R.RunData.TBDIS = dCov.TBDIS_s(:,i);
        R.Fit;
        FitResults{i} = R.FitResult;
        mNuSq(i) = R.FitResult.par(1);
    end
    %%
    MakeDir(savedir);
    save(savefileResult,'mNuSq','FitResults');
end



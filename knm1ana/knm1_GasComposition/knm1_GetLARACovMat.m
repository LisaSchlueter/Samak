% compute gas composition covariance matrix

nSamples = 5e3;
savedir  = [getenv('SamakPath'),'knm1ana/knm1_GasComposition/results/'];
savefileCovMat = sprintf('%sknm1_LARADataCovMat_%.0fsamples.mat',savedir,nSamples);
if exist(savefileCovMat,'file')
    dCov = importdata(savefileCovMat);
else
    
    % Init Model Object
    R = MultiRunAnalysis('RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1,...
        'minuitOpt','min ; minos',...
        'FSDFlag','SibilleFull',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF');
    R.exclDataStart = R.GetexclDataStart(40);  % 40eV range = 27 subruns
    R.InitModelObj_Norm_BKG;
    R.ModelObj.BKG_RateSec_i = 1e-05; % set background to essentially zero
    R.ModelObj.ComputeTBDDS;
    R.ModelObj.ComputeTBDIS;
    TBDIS_i = R.ModelObj.TBDIS;
    %% GET CORRELATION matrix
    savefileCorrMat = sprintf('%sknm1_LARAData.mat',savedir);
    d = importdata(savefileCorrMat);
    
    % input: trueness values from https://www.mdpi.com/1424-8220/20/17/4827
    WGTS_MolFrac_TT_RelErr = 9.7e-03;
    WGTS_MolFrac_HT_RelErr = 1.3e-01;
    WGTS_MolFrac_DT_RelErr = 9.2e-02;
    
    WGTS_MolFrac_TT_AbsErr =  WGTS_MolFrac_TT_RelErr.*R.RunData.WGTS_MolFrac_TT;
    WGTS_MolFrac_HT_AbsErr =  WGTS_MolFrac_HT_RelErr.*R.RunData.WGTS_MolFrac_HT;
    WGTS_MolFrac_DT_AbsErr =  WGTS_MolFrac_DT_RelErr.*R.RunData.WGTS_MolFrac_DT;
    
    AbsErr_Std = [WGTS_MolFrac_TT_AbsErr,WGTS_MolFrac_HT_AbsErr,WGTS_MolFrac_DT_AbsErr];
    CovMatLARA = AbsErr_Std.*d.CorrMat.*AbsErr_Std';
    
    % sampling
    WGTS_MolFrac_s =  mvnrnd([R.RunData.WGTS_MolFrac_TT,R.RunData.WGTS_MolFrac_HT,R.RunData.WGTS_MolFrac_DT],...
        CovMatLARA,nSamples);
    
    TBDIS_s = zeros(R.ModelObj.nqU,nSamples);
    for i=1:nSamples
        progressbar(i./nSamples);
        R.ModelObj.WGTS_MolFrac_TT = WGTS_MolFrac_s(i,1);
        R.ModelObj.WGTS_MolFrac_HT = WGTS_MolFrac_s(i,2);
        R.ModelObj.WGTS_MolFrac_DT = WGTS_MolFrac_s(i,3);
        
        R.ModelObj.ComputeTBDDS;
        R.ModelObj.ComputeTBDIS;
        TBDIS_s(:,i) = R.ModelObj.TBDIS;
    end
   
    covmat = cov(TBDIS_s');
    covmatfrac = covmat./(TBDIS_i.*TBDIS_i');
 
    
    save(savefileCovMat,'covmat','covmatfrac','TBDIS_s','WGTS_MolFrac_s','CovMatLARA','AbsErr_Std','TBDIS_i');
end
%%
%  GetFigure
%  imagesc(dCov.covmatfrac(13:end,13:end)-covmatfracshape(13:end,13:end));
%  c = colorbar;
%% convert to shape only
if exist(savefileCovMat,'file') && ~isfield(dCov,'covmatshape')
    dCov = importdata(savefileCovMat);
    NormCorr = mean(sum(dCov.TBDIS_s(13:end,:)))./sum(dCov.TBDIS_s(13:end,:)); % bring all spectra to the same number of counts (total)
    TBDIS_NormCorr = dCov.TBDIS_s.*NormCorr;
    
    covmatshape = cov(TBDIS_NormCorr');
    covmatfracshape = covmatshape./(dCov.TBDIS_i.*dCov.TBDIS_i');
    save(savefileCovMat,'covmatshape','covmatfracshape','-append');
else
    fprintf('not available \n');
     return
end



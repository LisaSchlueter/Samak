% knm2: parameter correlations
% based on fits of randomized twins

% load fit results
nFit = 2000;
RecomputeFlag = 'OFF';
range     = 40;
freePar   = 'mNu E0 Bkg Norm';
chi2      = 'chi2CMShape';
DataType  = 'Twin';
AnaFlag   = 'StackPixel';%'Ring';
RingMerge = 'None';
DopplerEffectFlag = 'FSD';
BKG_PtSlope = 3*1e-06;
TwinBias_BKG_PtSlope = 3*1e-06;
FSDFlag   = 'KNM2_0p1eV';
SysBudget    = 40;
NonPoissonScaleFactor = 1.112;
SigmaSq =  0.0124+0.0025;
savedir = [getenv('SamakPath'),'knm2ana/knm2_Correlations/results/'];
savefile = sprintf('%sknm2FitRandData_%s_%s_NP%.4g_%s_%.0feV_%s_%.0ffit.mat',...
    savedir,DataType,chi2,NonPoissonScaleFactor,strrep(freePar,' ',''),range,FSDFlag,nFit);

if strcmp(chi2,'chi2CMShape')
    savefile = strrep(savefile,chi2,sprintf('%s_SysBudget%.0f',chi2,SysBudget));
end


if exist(savefile,'file')
    d = importdata(savefile);
    d.FitPar(2,:) = d.FitPar(2,:)+d.A.ModelObj.Q_i; % endpoint
    d.FitPar(3,:) = (d.FitPar(3,:)+d.A.ModelObj.BKG_RateSec_i).*1e3; % background
     d.FitPar(4,:) = d.FitPar(4,:)+1;% signal normalization

    fprintf('load from file %s \n',savefile)
else
    fprintf(2,'file not found: %s \n',savefile)
    return
end

%%
corrplot_Uniform(d.FitPar(1:4,:))
 
[MeanCorrMat,StdCorrMat,CorrMat_i] = corrcoeff_err(d.FitPar(1:4,:));
%
 pltdir  = [getenv('SamakPath'),'knm2ana/knm2_Correlations/plots/'];
 MakeDir(pltdir);
 pltname = sprintf('%sknm2_Correlations.pdf',pltdir);
 %print(pltname,'-dpng','-r400');
 export_fig(pltname);
 fprintf('save plot to %s \n',pltname);
 
 

% knm1: parameter correlations
% based on fits of randomized twins


% load fit results
savedir = [getenv('SamakPath'),'knm1ana/knm1_FitRandData/results/'];
savefile = sprintf('%sknm1FitRandData_%s_%s_SysBudget24_NP%.4g_%s_%.0feV_%s_TwinmNuSq%.3geV2_Combi2000fits.mat',...
    savedir,'Twin','chi2CMShape',1.064,'mNuE0NormBkg',40,'Sibille0p5eV',-0.97);


if exist(savefile,'file')
    d = importdata(savefile);
    d.FitPar(2,:) = d.FitPar(2,:)+d.Real.ModelObj.Q_i; % endpoint
    d.FitPar(3,:) = (d.FitPar(3,:)+d.Real.ModelObj.BKG_RateSec_i).*1e3; % background
     d.FitPar(4,:) = d.FitPar(4,:)+1;% signal normalization

    fprintf('load from file %s \n',savefile)
else
    fprintf(2,'file not found: %s \n',savefile)
    return
end

%%
 corrplot_Uniform(d.FitPar(1:4,:))
 
 pltdir  = [getenv('SamakPath'),'knm1ana/knm1_FitRandData/plots/'];
 MakeDir(pltdir);
 pltname = sprintf('%sknm1_Correlations.png',pltdir);
 print(pltname,'-dpng','-r400');
% export_fig(pltname);
 fprintf('save plot to %s \n',pltname);
% knm1 original settings
% load single run fit results via object (uniform, fixed nu-mass)
% with NP factor
close all

savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits_NP.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

CommonArg = {'HideGaps','OFF','ShowRWPeriods','OFF','saveplot','ON','TimeStyle','date'};
% %R.PlotFitRunList(CommonArg{:},'Parameter','RhoD','YLim',1e17.*[1.076,1.15]);
% %R.PlotFitRunList('Parameter','T2');
%  R.PlotFitRunList('Parameter','B',CommonArg{:})
% R.PlotFitRunList('Parameter','pVal',CommonArg{:})
% FitResults = R.PlotFitRunList('Parameter','E0',CommonArg{:},'YLim',[18572.6 18575.2]);
% 
% %R.PlotFitRunList('HideGaps','OFF','ShowRWPeriods','OFF','YLim',[18572.6,18575.3],'Parameter','E0','saveplot','pdf');
% %R.PlotFitRunList('HideGaps','OFF','ShowRWPeriods','OFF','Parameter','B','saveplot','pdf');
% R.PlotFitRunList('Parameter','N',CommonArg{:},'DisplayStyle','Rel');%YLim',1+[-1,+1].*0.12
% 

%% evaluate stability
ndof = numel(FitResults.E0)-1;
chi2E0 = (wmean(FitResults.E0,1./FitResults.E0Err.^2)-FitResults.E0).^2./FitResults.E0Err.^2;
pE0 = 1-chi2cdf(sum(chi2E0),ndof);
[hE0, pE0r] =runstest(FitResults.E0); % runs test test for randomness, h = 0 -> accept null hypothesis, that x is in random order

chi2N = (wmean(FitResults.N,1./FitResults.NErr.^2)-FitResults.N).^2./FitResults.NErr.^2;
pN = 1-chi2cdf(sum(chi2N),ndof);
[hN, pNr] =runstest(FitResults.N); % runs test test for randomness


chi2B = (wmean(FitResults.B,1./FitResults.BErr.^2)-FitResults.B).^2./FitResults.BErr.^2;
pB = 1-chi2cdf(sum(chi2B),ndof);
[hB, pBr] =runstest(FitResults.B); % runs test test for randomness

%% p-values

[hp, ppr] =runstest(FitResults.pValue); % runs test test for randomness
uni_dist = makedist('Uniform');
[ksp,ksh]=kstest(FitResults.pValue,'CDF',uni_dist);
%mean(FitResults.E0Err)


cdf_s = cumsum(sort(FitResults.pValue))./sum(FitResults.pValue);



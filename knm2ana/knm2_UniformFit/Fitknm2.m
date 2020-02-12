Parameter='pVal';
knm2_UniformFit;
A.PlotFitRunList('Parameter',Parameter);
switch Parameter
    case 'E0'
        y=A.SingleRun_FitResults.chi2Stat.E0;
        yErr=A.SingleRun_FitResults.chi2Stat.E0Err;
    case 'B'
        y=A.SingleRun_FitResults.chi2Stat.B;
        yErr=A.SingleRun_FitResults.chi2Stat.BErr;
    case 'N'
        y=1+A.SingleRun_FitResults.chi2Stat.N;
        yErr=A.SingleRun_FitResults.chi2Stat.NErr;
    case 'pVal'
        y=A.SingleRun_FitResults.chi2Stat.pValue;
        yErr=A.SingleRun_FitResults.chi2Stat.pValueErr;
end
mean=wmean(y,1./yErr.^2);
s=std(y);
sErr=std(yErr);
wErr=ones(size(yErr));
Emean=wmean(yErr,wErr);
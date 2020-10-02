nGridSteps = 25;
chi2       = 'chi2CMShape';
freePar    = 'E0 Bkg Norm';
range      = 40;


KSNXGridSearch('nGridSteps',nGridSteps,...
    'chi2',chi2,...
    'freePar',freePar,...
    'range',range);
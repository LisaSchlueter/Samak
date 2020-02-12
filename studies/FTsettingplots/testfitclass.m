A = ref_FTMTD;
A.ComputeTBDDS();
A.ComputeTBDIS();
Data = {A.qU,A.TBDIS,A.TBDISE};

F = FITC('SO',A,'DATA',Data,'fitter','matlab','chi2name','chi2');

F.RESULTS

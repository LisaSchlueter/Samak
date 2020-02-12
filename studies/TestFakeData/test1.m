for i=1:1:10
    SPR.chi2 = 'chi2Stat';
    SPR.ModelObj.E0_bias=0;
    SPR.ModelObj.Bkg_Bias=0;
    SPR.ModelObj.Norm_Bias=0;
    SPR.ModelObj.mnuSq_Bias=0;
    SPR.Fit;
    disp(SPR.FitResult)
    
    Tpar(i,:) = SPR.FitResult.par;
    Terr(i,:) = SPR.FitResult.err;
    Tchi2(i)  = SPR.FitResult.chi2min;
    hold on
    SPR.PlotFitm2e0Contours
end
hold off

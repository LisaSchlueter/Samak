function X = fit_chi(m,s,R,eVr)
    % Return the X2 of a basic fit
    
    %R.ModelObj.mnuSq_i=-0.956847;
    R.ModelObj.mnu4Sq_i=m;
    R.ModelObj.sin2T4_i=s;
    R.exclDataStart = R.GetexclDataStart(eVr);
    R.RunData;
    R.Fit;
    
    X = R.FitResult.chi2min;
    
end
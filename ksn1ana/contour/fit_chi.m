function X = fit_chi(m,s,R,eVr)
    % Return the X2 of a basic fit
    
    % Initialisation
    R.i_Q=0;
    R.i_B=[];
    R.i_N=[];
    
%     R.ModelObj.mnuSq_i=-0.956847;                 % Neutrino mass
    R.ModelObj.mnu4Sq_i=m;                          % Sterile mass
    R.ModelObj.sin2T4_i=s;                          % Mixing angle
    
    R.exclDataStart = R.GetexclDataStart(eVr);
    R.RunData;
    R.Fit;
    
    X = R.FitResult.chi2min;
    
end
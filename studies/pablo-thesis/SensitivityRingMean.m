

RingSys = load('data/SensitivityRINGSYS');
RingStat = load('data/SensitivityRINGSTAT');
close all;
[par,err,chi2min,dof] = SysUncertainty(RingSys.sigma_mc',RingStat.sigma_mc',RingSys.mu_mc');
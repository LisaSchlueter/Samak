function [par,err,chi2min,dof] = SysUncertainty(TotUnc,StatUnc,parRun)
%addpath(genpath('../../../Samak2.0'));
% TotUnc, StatUnc, parRun

sysUncRun = sqrt(TotUnc.^2 - StatUnc.^2);
CMsys = sysUncRun*sysUncRun';
CMsys = 0.8*CMsys + diag(StatUnc.^2) + (0);


fprintf('----------BEGIN FIT MINUIT-------------- \n');
parMean = mean(parRun);
tmparg = sprintf(['set pri -10;migrad minos']);
Args = {parMean, {parRun, CMsys}, '-c', tmparg};
[par, err, chi2min, errmat] = fminuit('chi2meanSys',Args{:});

dof = length(parRun)-1;

results = {par, err, chi2min, errmat, dof};
end

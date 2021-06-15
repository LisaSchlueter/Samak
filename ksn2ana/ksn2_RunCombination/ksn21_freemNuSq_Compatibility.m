% compatibility test
% combine ksn1 and ksn2, nu-mass free
% plot with ksn1, ksn2, ksn1+2

% load mini combi file
savedir = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
MakeDir(savedir);
savefile = sprintf('%sksn21_Combi_freemNuSq_ReAna.mat',savedir);
load(savefile);

%%
Chi2Hat = d.chi2_ref - chi2ref_k1 -chi2ref_k2;
pHat = 1-chi2cdf(Chi2Hat,3);

fprintf(' Results --------------\n');
fprintf('KSN-1:  chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2  sin2t4=%.3f \n',chi2ref_k1,dof1,1-chi2cdf(chi2ref_k1,dof1),mNu4Sqbf_k1,sin2T4bf_k1);
fprintf('KSN-2:  chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2   sin2t4=%.3f \n',chi2ref_k2,dof2,1-chi2cdf(chi2ref_k2,dof2),mNu4Sqbf_k2,sin2T4bf_k2)
fprintf('KSN-12: chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2  sin2t4=%.3f \n',d.chi2_ref,dof12,1-chi2cdf(d.chi2_ref,dof12),mNu4Sqbf_k12,sin2T4bf_k12);
fprintf('KSN-12: chiHat^2=%.2f (%.0f dof) p = %.2f \n',Chi2Hat,3,pHat);

%fprintf('Compare to Null hypothesis DeltaChi^2 = %.2f (%.1f C.L.)\n',DeltaChi2_tot,100*chi2cdf(DeltaChi2_tot,2));


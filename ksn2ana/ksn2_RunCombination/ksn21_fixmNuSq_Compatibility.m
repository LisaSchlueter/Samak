% compatibility test
% combine ksn1 and ksn2, nu-mass fix
% plot with ksn1, ksn2, ksn1+2

DataType = 'Real';
% load mini combi file
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_RunCombination/results/'];
MakeDir(savedir)
savefile = sprintf('%sksn21_Combination_ReAna_%s.mat',savedir,DataType);
load(savefile);
%%
Chi2Hat = chi2_ref_12-chi2_ref_1-chi2_ref_2;
pHat = 1-chi2cdf(Chi2Hat,2);
fprintf(' Results --------------\n');
fprintf('KSN-1:  chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2  sin2t4=%.3f \n',chi2_ref_1,dof1,1-chi2cdf(chi2_ref_1,dof1),mNu4Sq_bf_1,sin2T4_bf_1);
fprintf('KSN-2:  chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2   sin2t4=%.3f \n',chi2_ref_2,dof2,1-chi2cdf(chi2_ref_2,dof2),mNu4Sq_bf_2,sin2T4_bf_2)
fprintf('KSN-12: chi^2=%.2f (%.0f dof) p = %.2f ,  m4^2 = %.1f eV^2  sin2t4=%.3f \n',chi2_ref_12,dof12,1-chi2cdf(chi2_ref_12,dof12),mNu4Sq_bf_12,sin2T4_bf_12);
fprintf('KSN-12: chiHat^2=%.2f (%.0f dof) p = %.2f \n',Chi2Hat,2,pHat);

% fprintf('Compare to Null hypothesis DeltaChi^2 = %.2f (%.1f C.L.)\n',...
%     chi2_null_12-chi2_ref_12,100*chi2cdf(chi2_null_12-chi2_ref_12,2));
 
% dof for combi map: fixed m^2
% KSN-1:  dof = 27(nqU)-5(sin2T4,m4Sq,E0,Norm,Bkg) = 22;
% KSN-2:  dof = 28(nqU)-5(sin2T4,m4Sq,E0,Norm,Bkg) = 23;
% KSN-12: dof = 55(nqU)-5(sin2T4,m4Sq,E01,Norm1,Bkg1,E02,Norm2,Bkg2) = 47;
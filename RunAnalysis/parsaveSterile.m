function parsaveSterile(fname,FitResults_Null,FitResults_NoSterile,chi2_ref)
save(fname,'FitResults_Null','FitResults_NoSterile','chi2_ref','-append');
fprintf('Save fit result to file %s \n',fname)
end
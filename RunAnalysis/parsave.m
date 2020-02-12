function parsave(fname, FitResult,Run,FitListArg)
save(fname,'FitResult','Run','FitListArg');
 fprintf('Save fit result to file %s \n',fname)
end
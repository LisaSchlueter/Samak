% test rebinned fsd

% rebinned fsd: excited states rebined, ground state not rebinned
RunList = 'KNM1';
A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','FSDFlag','Sibille0p5eV','exclDataStart',14);
tic; A.Fit; toc

% default (finder binned) fsd
B = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','FSDFlag','Sibille','exclDataStart',14);
tic; B.Fit; toc

%%
fprintf(' ------------------------------------ \n');
fprintf('difference nu-mass = %.3g \n',A.FitResult.par(1)-B.FitResult.par(1));
fprintf('difference nu-mass err = %.3g \n',A.FitResult.err(1)-B.FitResult.err(1));
fprintf(' ------------------------------------ \n');
fprintf('difference E0 = %.3g \n',A.FitResult.par(2)-B.FitResult.par(2));
fprintf('difference E0 err = %.3g \n',A.FitResult.err(2)-B.FitResult.err(2));
fprintf(' ------------------------------------ \n');
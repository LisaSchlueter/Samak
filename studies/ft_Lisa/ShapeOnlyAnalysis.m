addpath(genpath('../../../Samak2.0'));
RunNrs = [40531, 40538,40539,40540,40541,40542,40543];

% How we always do it....
Ana = RunAnalysis('RunNr',RunNrs,'AnaFlag','StackPixel','chi2','chi2CM');
Ana.Fit
ResultsCM = Ana.FitResult;

AnaNormfix = RunAnalysis('RunNr',RunNrs,'AnaFlag','StackPixel','chi2','chi2CM');
AnaNormfix.fixPar ='4';
AnaNormfix.Fit;
ResultsNormFix = AnaNormfix.FitResult;

AnaShape = RunAnalysis('RunNr',RunNrs,'AnaFlag','StackPixel','chi2','chi2CMShape');
AnaShape.Fit;
ResultsCMShape = AnaShape.FitResult;



%% plots
meanE0 = mean([ResultsCM.par(2), ResultsCMShape.par(2), ResultsNormFix.par(2)]);
meanB = mean([ResultsCM.par(3), ResultsCMShape.par(3), ResultsNormFix.par(3)]);
meanN = mean([ResultsCM.par(4), ResultsCMShape.par(4), ResultsNormFix.par(4)]);
meanchi2 = mean([ResultsCM.chi2min, ResultsCMShape.chi2min, ResultsNormFix.chi2min]);

subplot(2,2,1); %E0
bar(1:3,[ResultsCM.par(1), ResultsNormFix.par(1), ResultsCMShape.par(1)]);
title('E0 (eV)');

subplot(2,2,2); %Background
bar(1:3,[ResultsCM.par(3), ResultsNormFix.par(3), ResultsCMShape.par(3)]);
title('BKG (cps)');

subplot(2,2,3); %Norm
bar(1:3,[ResultsCM.par(4), ResultsNormFix.par(4), ResultsCMShape.par(4)]);
title('N');

subplot(2,2,4); %chi2
bar(1:3,[ResultsCM.chi2min, ResultsNormFix.chi2min, ResultsCMShape.chi2min]);
title('\chi^2');






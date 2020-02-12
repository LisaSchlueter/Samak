%------------------------------------------------------------------------%
% Script to Compute Endpoint Results for Samak Summary
% Lisa Schlueter
%% ------------------------------------------------------------------------%

RunList = [40531,40538:40543,40603,40604,40610:40613,40667:40693];
chi2Flag = 'chi2Stat';
exclDataStart = 9;
myResults = FitRunList2('RunList',RunList,'ringCutFlag','ex2','AnaFlag','StackPixel',...
           'fixPar','1 5 6','chi2',chi2Flag,'exclDataStart',exclDataStart);
close all

%% Write Fit Results to File
cd plots/tmp
filename = sprintf('FitResults_%s_exclDataStart%u_%uRuns_%u_%u.txt',chi2Flag,exclDataStart,numel(RunList),RunList(1),RunList(end));
unix(['touch ',filename]);
Chi2 = [myResults.chi2min; myResults.dof; myResults.pValue];
Endpoint = [myResults.E0+18575; myResults.E0Err];
Background = [myResults.B; myResults.BErr];
Normalization = [myResults.N+1; myResults.NErr];

fileID = fopen(filename,'w');
mytitle = sprintf('Fit to %u Runs: %u-%u \n',numel(RunList),RunList(1),RunList(end));
fprintf(fileID,mytitle);
fprintf(fileID,'mean E0 =%.3f \n',mean(myResults.E0+18575));
fprintf(fileID,'mean E0err = %.3f \n', mean(myResults.E0Err));
fprintf(fileID,'err of mean (E0) = %.3f \n', sqrt(var(myResults.E0)));
fprintf(fileID,'chi2/dof \n');
fprintf(fileID,'%.2f / %.0f \n',[myResults.chi2min; myResults.dof]);
fprintf(fileID,'p-value \n');
fprintf(fileID,'%.3f \n',myResults.pValue);
fprintf(fileID,'Endpoint (eV) \n');
fprintf(fileID,'%.3f +/- %.3f \n',Endpoint);
fprintf(fileID,'Backgroud (mps) \n');
fprintf(fileID,'%.3f +/- %.3f \n',Background);
fprintf(fileID,'Normalization \n');
fprintf(fileID,'%.4f +/- %.4f \n',Normalization);
fclose(fileID);
cd ../..

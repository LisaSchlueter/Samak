addpath(genpath('../../../Samak2.0'));
RunList = [40531,40538:40543,40603,40604,40610:40613,40667:40693];
MRA = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','chi2','chi2CMShape',...
'exclDataStart',1);
MRA.RhoDScan('saveplot','ON');
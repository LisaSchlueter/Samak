addpath(genpath('../../../Samak2.0'));
%R = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Aseev');
%R.ComputeCM('nTrials',1000,'Stack','ON','SysEffect',struct('RF_BF','ON','RF_EL','ON'));
R = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Abdurashitov');
R.ComputeCM('nTrials',1000,'Stack','OFF');
R = MultiRunAnalysis('RunList','StackCD100all','chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Aseev');
R.ComputeCM('nTrials',1000,'Stack','OFF');

R = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Abdurashitov');
R.ComputeCM('nTrials',1000,'Stack','OFF');
R = MultiRunAnalysis('RunList','StackCD100_3hours','chi2','chi2Stat','DataEffCor','RunSummary','ELossFlag','Aseev');
R.ComputeCM('nTrials',1000,'Stack','OFF');



%R = MultiRunAnalysis('RunList','StackCD100up','chi2','chi2Stat','DataEffCor','RunSummary');
% R.ComputeCM('nTrials',5000,'Stack','ON');
% % R = MultiRunAnalysis('RunList','StackCD100down','chi2','chi2Stat','DataEffCor','RunSummary');
% % R.ComputeCM('nTrials',5000,'Stack','ON');
% 


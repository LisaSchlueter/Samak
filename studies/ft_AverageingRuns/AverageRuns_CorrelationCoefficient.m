%-----------------------------------------------------------------------------------------------------
% Function to estimate correlation coefficient of systematics
% Compute Mean Endpoint (27 Runs)
% Compute Standard Error of Mean
% Repeat as a function of correlation coefficient
% 
% L.Schlueter (09/2018)
%-----------------------------------------------------------------------------------------------------
addpath(genpath('../../../Samak2.0'));
%TD = 'StackCD100_3hours';
TD = 'StackCD100all';
exclDataStart = 7;
CorrCoeff = (0:0.1:1);
%Init
MeanE0       = zeros(numel(CorrCoeff),1);
ErrMeanE0    = zeros(numel(CorrCoeff),1);
chi2min      = zeros(numel(CorrCoeff),1);
dof          = zeros(numel(CorrCoeff),1);
SysCM        = cell(numel(CorrCoeff),1);
ResultsCMall = cell(numel(CorrCoeff),1);

MRA = MultiRunAnalysis('RunList',TD,'fixPar','1 5 6','exclDataStart',exclDataStart,...
    'DataEffCorr','ROI+PileUp','ringCutFlag','ex2','chi2','chi2CM');
progressbar('AverageRuns correlation coefficient script');
for i=1:numel(CorrCoeff)
    progressbar(i/numel(CorrCoeff));
    [MeanE0(i), ErrMeanE0(i), chi2min(i), dof(i), SysCM{i}, ResultsCMall{i}] = ...
        MRA.FitRunList_AveragePar('CorrCoeff',CorrCoeff(i),'saveplot','ON','ReFit','ON');
end
belowE0 = sprintf('%.0feVbelowE0',18575-MRA.ModelObj.qU(exclDataStart));
savename = sprintf('./results/AverageRuns_CorrelationCoefficient_RoiPileUp_%sex2_%s.mat',TD,belowE0);
save(savename,'MeanE0','ErrMeanE0','CorrCoeff','chi2min','dof','SysCM','ResultsCMall');
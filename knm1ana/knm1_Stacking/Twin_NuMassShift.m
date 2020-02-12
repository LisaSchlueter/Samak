% Short Script to fit twins switch and without systematics
% investigate neutrino mass shift and endpint shift
RunList = 'KNM1';
chi2name = 'chi2Stat';

exclDataStartAll = [17,14,2];
nRange = numel(exclDataStartAll);
mNu    = zeros(nRange,1);
mNuErr = zeros(nRange,1);
E0 = zeros(nRange,1);
E0Err= zeros(nRange,1);

if ~exist('T','var')
    T = MultiRunAnalysis('RunList',RunList,'DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11','chi2',chi2name);
end

 T.chi2 = chi2name;
if ~strcmp(chi2name,'chi2Stat')
    T.ComputeCM('SysEffects',struct('Stack','ON'),'BkgCM','OFF');
end

for i=1:nRange
T.exclDataStart=exclDataStartAll(i);
T.Fit;
mNu(i) = T.FitResult.par(1); mNuErr(i) = T.FitResult.err(1);
E0(i)  = T.FitResult.par(2); E0Err(i)  = T.FitResult.err(2);
end


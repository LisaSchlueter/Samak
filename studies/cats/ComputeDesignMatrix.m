function DesignMatrix=ComputeDesignMatrix(f,FitParBest,varargin)
%
% Compute Design Matrix https://en.wikipedia.org/wiki/Design_matrix
% Input: 
%  f: function
%  FitPar: fit parameters
%  FitParList
%
% Output:
%  DesignMatrix
%
% T. Lasserre, June 2018
%

% Parser
p = inputParser;
p.addParameter('epsilon',1e-6,@(x)isfloat(x));
p.addParameter('FitParList',1:numel(FitParBest));
p.parse(varargin{:});
epsilon       = p.Results.epsilon;
FitParList    = p.Results.FitParList;

% Init
nFitPar      = length(FitParBest);  
f0           = feval(f,FitParBest);
nBins        = length(f0);  
DesignMatrix = zeros(nBins,nFitPar);

% Compute Design Matrix
for i=FitParList
    Variation = [zeros(i-1,1); epsilon; zeros(nFitPar-i,1)];
    DesignMatrix(:,i) = ( feval(f,FitParBest+Variation) - feval(f,FitParBest-Variation)) / (2*epsilon);
end
    




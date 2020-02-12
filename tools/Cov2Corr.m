%% Converts covariance to variance and correlation coefficient
function [ExpVar, ExpCorrC] = Cov2Corr(ExpCovariance)

%   Cov2Corr Converts covariance to variance and correlation coefficient.
%   Computes the volatilities of N random processes and the degree of
%   correlation between the processes.  
%
%   [ExpVar, ExpCorrC] = cov2corr(ExpCovariance)
%
%   Input:
%     ExpCovariance: N by N covariance matrix, e.g. from COV or EWSTATS
%
%   Outputs:
%     ExpVar     : 1 by N vector with the variances of each process
%
%     ExpCorrC     : N by N matrix of correlation coefficients.  The
%     entries of ExpCorrC range from 1 (completely correlated) to -1
%     (completely anti-correlated).  A value of 0 in the (i,j) entry
%     indicates that the i'th and j'th processes are uncorrelated.
% 
%   Expvar(i) = ( ExpCovariance(i,i) );
%   ExpCorrC(i,j) = ExpCovariance(i,j)/sqrt( ExpVar(i)*ExpVar(j) );
% 
%   See also EWSTATS, COV, CORRCOEF, VAR, CORR2COV.

%-----------------------------------------------------------------
% Argument checking
% ExpCovariance   [N by N]  with diag(ExpCovariance)>=0
% N      [scalar]
%-----------------------------------------------------------------
if nargin < 1,
  error('finance:Cov2Corr:missingInput','Enter a covariance matrix.')
end

if size(ExpCovariance,1) ~= size(ExpCovariance, 2)
  error('finance:Cov2Corr:invalidCovMatrixSize','Covariance matrix must be square')
else
  N = size(ExpCovariance, 1);
end

if any( diag(ExpCovariance) < 0 )
  error('finance:Cov2Corr:invalidCovMatrixSymmetry','Covariance matrix must be symmetric with non-negative diagonal')
end

%-----------------------------------------------------------------
% Simple correlation is ExpCovariance./( ExpSigma'*ExpSigma )
% ExpSigma [1 by N]
% ExpCorrC [N by N]
%-----------------------------------------------------------------
ExpVar = (diag(ExpCovariance))';

% start with default correlation of identity for degenerate processes
ExpCorrC = eye(N);

% find processes which are not degenerate
IndPos = (ExpVar > 0);

% Compute correlation only for non-degenerate processes
ExpCorrC(IndPos,IndPos) = ExpCovariance(IndPos,IndPos) ./ ... 
    sqrt(ExpVar(IndPos)' * ExpVar(IndPos));

%-----------------------------------------------------------------
% end of function Cov2Corr
%-----------------------------------------------------------------

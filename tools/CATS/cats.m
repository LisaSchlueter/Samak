function stats = cats(y,A,varargin)
%CATS Chisquare Analysis Tool with Systematics
%
% stats = cats(y,A,options)
% - y is a [n,1] vector of data to be fitted
% - A is a [n,p] matrix (the design matrix) of the linear model with p
%   columns for the p parameters of the model.
% - stats is the output structure containing the fit outputs. Different
%   level of outputs and chi square treatment available through the options
%   provided to the function.
%
% 'options' are specified as follows:
%
% - 'Vinv',Vinv is a [n,n] matrix: this is the inverse covariance matrix of the data
% - 'S', S  is a [n,k] matrix (the systematic design matrix) of the linear
%   model with k columns for the k parameters of the model.
% - 'x0',x0 is a [k,1] vector of a priori values on known systematic parameters
% - 'Winv',Winv is a [k,k] matrix: this is the inverse covariance matrix of the known
%   systematic parameters.
% - 'L',L is a [c,p+k] matrix: this is the linear constraint matrix on the k
%   parameters.
% - 'r',r is a [c,1] vector of the c constraint on the p+k parameters.
%
% y and A are compulsory arguments. All the other arguments are optional.
% If not specified, default value for Vinv is the identity matrix (no a
% priori knowledge on the variance-covariance of the data). By default if
% not specified, S, x0, r and Winv, L are empty vectors and matrices.
%
% The function returns a structured object called here stats whose
% different fields bring information on the fit, its quality, the model
% parameters at minimum of the chi square and the associated uncertainties
% on parameters and their correlations. Moreover, the function returns also
% in stats structure a lot of information about the fit through diganosis
% variables.
%
% Diagnostic levels:
% - Level 0: Model coherence, reliability (model colinearities,...)
% - Level 1: Data/Model global coherence, agreement (chi2, residuals,...)
% - Level 2: Influencial/Atypical points detection (leverages, Cook's D
%            Covariance ratios,...)
% - Level 3: Plans to add for future goodness of fit tests and
%            compatibility tests among provided inputs: parameter goodness
%            of fit test, Fisher tests,...
% - Level 4: Pattern identification, beyond the chi2 (run test,
%            Durbin-Watson,...) residual test??? possible? - Check
%            Neyman-Barton smooth test test for H1, order k alternative as
%            with Rayner-Best tests
% - Level 5: DSD = Data Set Diagonalization: Pumplin et al.
% - Level 6: sPlot = Le Diberder et al.
%
% stats = CATS(y,A,...,'RenormalizeErrors',true) changes the renormalize
%  error option of the fitted parameter uncertainties. This is a global option
%  set up by default to false. If it is set up to true as in the above statement
%  the uncertainties and covariances on fitted parameters are multiplied by sqrt(mse)
%  (i.e. sqrt(chi2min/ndof)).
%
% stats = CATS(y,A,...,'Names',{'...'}) helps identifying the data input
% sources, as data points, systematic pulls through names given as a
% structure of strings. The number of elements in this structure should be
% the same as the number of bins plus the number of systematic parameters.
%
% stats = CATS(y,A,...,'whichstats',{'...'}) creates an output
%   structure stats containing the statistics listed in whichstats.
%   whichstas can be a single string such as 'leverage' or a cell array of
%   strings such as {'leverage' 'standres' 'studres'}.  By default, 
%   CATS returns all statistics. Valid statistic strings are:
%
%      Name          Meaning
%      'Q'           Q from the QR Decomposition of the design matrix
%      'R'           R from the QR Decomposition of the design matrix
%      'xhat'        Regression coefficients
%      'covx'        Covariance of regression coefficients
%      'yhat'        Fitted values of the response data
%      'r'           Residuals
%      'mse'         Mean squared error
%      'rsquare'     R-square statistic
%      'adjrsquare'  Adjusted R-square statistic
%      'leverage'    Leverage
%      'hatmat'      Hat (projection) matrix
%      's2_i'        Delete-1 variance
%      'beta_i'      Delete-1 coefficients
%      'standres'    Standardized residuals
%      'studres'     Studentized residuals
%      'dfxs'        Scaled change in regression coefficients
%      'dffit'       Change in fitted values
%      'dffits'      Scaled change in fitted values
%      'covratio'    Change in covariance
%      'cookd'       Cook's distance
%      'tstat'       t statistics for coefficients
%      'fstat'       F statistic
%      'dwstat'      Durbin Watson statistic
%                    Heteroscedasticity consistent covariance matrix tests
%      'hc0'         see Davidson & McKinnon
%      'hc1'
%      'hc2'
%      'hc3'
%      'hc4'
%      'hac'
%      'empty'
%      'rankdef'
%      'sqrtVinv'    square root of Vinv matrix
%      'sqrtWinv'    square root of Winv matrix
%      'corx'        Correlation matrix of coefficients
%      'pi'          Collinearity diagnostics: PI matrix
%      'heta'        Collinearity diagnostics: heta numbers, kth condition index
%      'all'         Create all of the above statistics

%   References:
%   Belsley, D.A., E. Kuh, and R.E. Welsch (1980), Regression
%      Diagnostics, New York: Wiley.
%   Cook, R.D., and S. Weisberg (1982), Residuals and Influence
%      in Regression, New York: Wiley.

%% CATS Todo List - August 2011
% Status: Level 1 and 2 implemented. Level 0 is ongoing developments. Level
% 3 is to be done - August 2011. 
%
% 1/ add the confidence/prediction intervals on the parameters/model
% (functionnal and observational prediction intervals)
% 
% 2/ Level 0 diagnostics (colinearity)
%   = Column deletion, VIF (Variance inflation factor, condition number)
%
% 3/ Level 3 diagnostics (beyond chi2 tests)
%   Multiple row deletion instead of single row deletion. Allow to group
%   information from multiple bins. Require a strategy on the group
%   selection, otherwise this leads to complex combinatorial problems and
%   leads to 2step (or even n-steps) chi square analyses for principal
%   component analyses => requires heavier developments. First
%   implementation steps could be to have groups of bin info provided by
%   the user...
%

% Work based on:
%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.5 $  $Date: 2011/12/27 16:17:02 $

% This code:
%   Copyright 2005-2014 Guillaume MENTION, CEA Saclay
%   $Revision: 1.32$  $Date: 2014/04/20$

%% Define useful sizes
nobs = size(y,1); % # of data points
npar = size(A,2); % # of physical parameters

%% Prepare the validity checking of the input arguments

% Create a parser
P = inputParser;
% Compulsory arguments
P.addRequired('y',@(x)isfloat(x) && ismatrix(x) && size(x,2)==1); % Data
P.addRequired('A',@(x)isfloat(x)); % Linear model
% Optional arguments
% Data
P.addParameter('Vinv',eye(nobs),@(x)isfloat(x) && size(x,1)==nobs && size(x,2)==nobs); % uncertainties on data
% Systematics
P.addParameter('S',[],@(x)isfloat(x) && ismatrix(x)); % Linear model matrix of systematics
P.addParameter('Winv',[],@(x)isfloat(x) && ismatrix(x) && size(x,1)==size(x,2)); % % prior uncertainties on systematics
P.addParameter('x0',[],@(x)isfloat(x) && ismatrix(x) && size(x,2)==1); % prior biases on systematics
% Constraints
P.addParameter('L',[],@(x)isfloat(x) && ismatrix(x)); % Linear constraints on model parameters
P.addParameter('r',[],@(x)isfloat(x) && ismatrix(x) && size(x,2)==1); % constraining values on parameters
P.addParameter('whichstats','all',@(x)ischar(x)||iscell(x));
% Parameter names (sys,par)
P.addParameter('RenormalizeErrors',false,@islogical);
P.addParameter('Names',{},@(x)all(cellfun(@ischar,x)));
%% Parse the input arguments
P.parse(y,A,varargin{:});

if ~isempty(A)&&~(ismatrix(A) && size(A,1)==nobs)
    error('CATS:core','Design matrix A has to be a nobs * npar matrix.');
end


%% Store the results
Vinv = P.Results.Vinv;

% Systematics argument checking
x0 = P.Results.x0;
Winv = P.Results.Winv;
S = P.Results.S;

%% Perform the last checkings
if isempty(S)&&(~isempty(Winv)||~isempty(x0))
    warning('CATS:core','Systematics not included since S matrix is not specified.');
    Winv = [];
    x0 = [];
    %nsys = 0;
else
    if ~isempty(S)&&~isempty(Winv)&&isempty(x0)
        if (size(Winv,1)~=size(S,2))
            error('CATS:core','The number of colums in S has to be the same as the number of columns and rows in Winv!');
        else
            nsys = size(Winv,1);
        end
        x0 = zeros(nsys,1);
    end
    if ~isempty(S)&&isempty(Winv)&&~isempty(x0)
        if (size(x0,1)~=size(S,2))
            error('CATS:core','The number of colums in S has to be the same as the number of rows in x0!');
        else
            nsys = size(x0,1);
        end
        Winv = eye(nsys);
    end
    if ~isempty(S)&&isempty(Winv)&&isempty(x0)
        warning('CATS:core','You did not provide any knowledge about the systematics, thus they are treated as physical paremeters.');
        A = [S A];
        S = [];
        npar = size(A,2);
        %nsys = 0;
    end
end

nsys = size(x0,1);

% Constraints argument checking
L = P.Results.L;
r = P.Results.r;

if isempty(L)&&~isempty(r)
    ncst = size(r,1);
    L=eye(ncst,nsys+npar);
end
if ~isempty(L)&&isempty(r)
    if size(L,2)~=nsys+npar
        error('CATS:core','You have to provide the constraints on all the parameter and systematics at the same time.');
    else
        ncst = size(L,1);
        r=zeros(ncst,1);
    end
end

ncst = size(r,1);

if (numel((P.Results.Names))~=(nobs+nsys+npar)) && ~isempty(P.Results.Names)
    error('CATS:core','You did not provide not the right amount of names for systematics and parameters');
else
    names  = P.Results.Names;
end


%% Remove missing values, if any

wasnan = isnan(y);
if npar
    wasnan = wasnan | any(isnan(A),2);
end
if nsys
    wasnan = wasnan | any(isnan(S),2);
end
havenans = any(wasnan);
if havenans
    y(wasnan) = [];
    if npar
        A(wasnan,:) = [];
    end
    if nsys
        S(wasnan,:) = [];
    end
    Vinv(wasnan,:) = [];
    Vinv(:,wasnan) = [];
    nobs = length(y);
end

%% Define basic and inputs values
stats.y = y;
stats.A = A;
stats.S = S;
stats.Vinv = Vinv;
stats.Winv = Winv;
stats.x0 = x0;
stats.L = L;
stats.c = r;

stats.names = names;

stats.nobs = nobs;
stats.nsys = nsys;
stats.npar = npar;
stats.ncst = ncst;

%% Reduce the LS problem to OLS
% sqrtVinv = chol(Vinv);
% Be careful eig changes row order!!!! Reduces residual indexes do not
% correspond to the intial ones... 07/05/2013
% Got the solution on April 2014!!!!! instead of sqrt(D)*P' use
% P*sqrt(D)*P' !!!!! Of course !!!!!!!!! any matrix function:
% f(M) = P*diag(f(diag(D)))*P^-1, P^-1 = P' for symmetric matrices !!!
[Pv,Dv] = eig(Vinv); % Be careful with definition: Vinv = Pv*Dv*Pv';
sqrtVinv = Pv*sqrt(Dv)*Pv'; % very important to put the Pv on the left to go
% back to the initial basis!!!!!! - GM.- April 2014 comment!!!!

% sqrtWinv = chol(Winv);
[Pw,Dw] = eig(Winv);
sqrtWinv = Pw*sqrt(Dw)*Pw';  % very important to put the Pw on the left to go
% back to the initial basis!!!!!! GM.- April 2014 comment!!!!

%% Generalized linear model including data correlations and systematics
%  and possible linear constraints on parameters

% Defining the generalized design matrix and the generalized data vector
if nsys
    if npar
        Ah = [sqrtVinv*S sqrtVinv*A; sqrtWinv zeros(nsys,npar)];
    else
        Ah = [sqrtVinv*S; sqrtWinv];
    end
else
    Ah = sqrtVinv*A;
end
yh = [sqrtVinv*y; sqrtWinv*x0];

if ncst
    Lp = pinv(L);
    yh = yh-Ah*Lp*r;
    Ah = Ah*(eye(npar+nsys)-Lp*L);
    % Equivalent other way to compute the restricted solution
    %N=null(L);
    %Ah = Ah*N*pinv(N'*N)*N';
    % Frederick James solution, but problem with the dimensions not solved
    % yet
    %Ah = [Ah pinv(Ah')*L';L zeros(length(r))];
    %yh = [yh; r];
end

% Accessible var names in the output structure
% add up 'corx' but be careful with indexes
varnames = {'Q','R','xhat','covx','yhat','r','mse','rsquare','adjrsquare', ...
    'leverage','hatmat','s2_i','x_i','standres','studres', ...
    'dfxs','dffit','dffits','covratio','cookd','tstat','fstat','dwstat',...
    'hc0','hc1','hc2','hc3','hc4','hac','empty','rankdef','sqrtVinv',...
    'sqrtWinv','corx','heta','pi'};

whichstats = P.Results.whichstats;

if ~iscell(whichstats), whichstats = {whichstats}; end
  if ~iscellstr(whichstats)
     error('CATS:BadStats',...
           'whichstas argument must be one or more statistic names.');
  end
  idx = false(length(varnames),1);
  for j=1:length(whichstats)
     snj = whichstats{j};
     if isequal(snj,'all')
        idx(:) = true;
        break;
     else
        k = find(strcmp(snj,varnames));
        if isempty(k)
           error('CATS:BadStats',...
                 'Invalid statistic name ''%s''.',snj);
        else
           idx(k) = true;
        end
     end
  end
  if ~any(idx), return, end
%idx = true(length(varnames),1);

%% Solve the Chi2

[Q,R] = qr(Ah,0);
if rank(Ah)<npar+nsys
    xhat = pinv(R)*Q'*yh;
else
    xhat = R\(Q'*yh);
end
if ncst
    xhat = xhat + Lp*r;
end

yrhat = Ah*xhat;

if nsys
    % yhat = [inv(sqrtVinv) zeros(nobs,nsys+npar);...
    %         zeros(nsys+npar,nobs) inv(sqrtWinv)]*yrhat;
    yhat =  [S A;eye(nsys) zeros(nsys,npar)]*xhat;
else
    yhat = sqrtVinv\yrhat;
end

%% Diagnoses computations

% OLS-form residuals with "decorrelated/deweighted" data
% OLS: Ordinary Least Squares
residuals = yh-yrhat;
% degrees of freedom of the chi2
dfe = nobs-npar+ncst;
stats.ndof = dfe;
% Estimation of variance of data implies the degrees of freedom of the total sum of squares
dft = nobs+nsys-1;
% Mean of the data
ybar = mean(yh);
% Sum of Squared Errors, Estimated Sum of Squares
sse = norm(residuals)^2;
% Refression Sum of Squares, or Fitted Sum of Squares
ssr = norm(yrhat - ybar)^2;
% Total Sum of Squares: sst = ssr + sse
sst = norm(yh - ybar)^2;
% Chi2min/ndof, the minimum chi2 per degree of freedom
mse = sse./dfe;
% The leverages, diagonal elements of the hat projection matrix
h = sum(abs(Q).^2,2);
% Leave 1 variances: computation of variance if the corresponding data
% point is removed
s_sqr_i = (dfe*mse - abs(residuals).^2./(1-h))./(dfe-1);
% Internally Studentized Residuals
e_i = residuals./sqrt(s_sqr_i.*(1-h));
e_i = e_i.*(abs(1-h)>eps); % avoid problems with full leverage points
% Variance/Covariance matrix of the parameters from the fit
if ncst
    ri = pinv(R);
else
    ri = R\eye(nsys+npar);
end
xtxi = ri*ri';
if (abs(mse)<1e-3)||isinf(mse)||~P.Results.RenormalizeErrors
    % sse
    covx = xtxi;
else
    covx = xtxi*mse;
end


%% TODO: add the confidence/prediction intervals on the parameters/model
% (functionnal and observational prediction intervals)
%
% Parameter confidence intervals
% if (nobs-npar ~= 0)
%     rmse = normr/sqrt(nu);    % Root mean square error.
%     tval = tinv((1-alpha/2),nu);
% else
%     rmse = NaN;
%     tval = 0;
% end
% se = zeros(ncolX,1);
% se(perm,:) = rmse*sqrt(sum(abs(rI).^2,2));
% xhat_int = [xhat-tval*se, xhat+tval*se];

%%

% Do one preliminary calculation
if idx(13) || idx(16)
    % Delete 1 coefficients. X_I
    stde = residuals./(1-h);
    stde = stde(:,ones(nsys+npar,1));
    x_i = xhat(:,ones(nobs+nsys,1)) - ri*(Q.*stde)';
end

% Store each requested statistic into the structure
% stats.source = 'regstats';
if idx(1)  % Q from the QR decomposition of the Ah matrix.
    stats.(varnames{1}) = Q;
end
if idx(2)  % R from the QR decomposition of the Ah matrix.
    stats.(varnames{2}) = R;
end
if idx(3)  % Coefficients.
    stats.(varnames{3}) = xhat;
end
if idx(4)   % Covariance of the parameters.
    stats.(varnames{4}) = covx;
end
if idx(5)  % Fitted values.
    tmp = yhat;
    if havenans, tmp = fixrows(tmp, wasnan); end
    stats.(varnames{5}) = tmp;
end
if idx(6)  % Residuals.
    tmp = [sqrtVinv zeros(size(sqrtVinv,1),size(sqrtWinv,2));
        zeros(size(sqrtWinv,1),size(sqrtVinv,2)) sqrtWinv]\residuals;
    if havenans, tmp = fixrows(tmp, wasnan); end
    stats.(varnames{6}) = tmp;
end
if idx(7)  % Mean squared error.
    stats.(varnames{7}) = mse;
end
if idx(8)  % R-square.
    % There are several ways to compute R^2, all equivalent for a linear
    % model where Ah includes a constant term, but not equivalent
    % otherwise.  R^2 can be negative for models without an intercept.
    % This indicates that the model is inappropriate.
    stats.(varnames{8}) = 1 - sse ./ sst;
end
if idx(9)  % Adjusted R-square.
    stats.(varnames{9}) = 1 - (sse./sst)*(dft./dfe);
end
if idx(10)  % Leverage.
    tmp = h;
    if havenans, tmp = fixrows(tmp, wasnan); end
    stats.(varnames{10}) = tmp;
end
if idx(11)  % Hat Matrix.
    H = Q*Q';
    if havenans
        tmp = zeros(length(wasnan));
        tmp(~wasnan,~wasnan) = H;
        tmp(wasnan,wasnan) = diag(NaN(sum(wasnan),1));
        H = tmp;
    end
    stats.(varnames{11}) = H;
end
if idx(12) % Delete 1 variance (mse). S_I
    tmp = s_sqr_i;
    if havenans, tmp = fixrows(tmp, wasnan); end
    stats.(varnames{12}) = tmp;
end
if idx(13) % Delete 1 coefficients. X_I
    tmp = x_i;
    if havenans
        % Estimates would be same if missing observations left out
        tmp = zeros(npar,length(wasnan));
        tmp(:,~wasnan) = x_i;
        tmp(:,wasnan) = xhat(:,ones(sum(wasnan),1));
    end
    stats.(varnames{13}) = tmp;
end
if idx(14) % Standardized residuals.
    standr = residuals./(sqrt(mse*(1-h)));
    if havenans, standr = fixrows(standr, wasnan); end
    standr(isnan(standr)) = 0;
    standr = standr.*(abs(1-h)>eps);
    stats.(varnames{14}) = standr;
end
if idx(15) % Studentized residuals.
    studr = e_i;
    if havenans, studr = fixrows(studr, wasnan); end
    studr(isnan(studr)) = 0;
    stats.(varnames{15}) = studr;
end
if idx(16) % Scaled change in xhat. DFXHAT
    x = xhat(:,ones(nobs+nsys,1));
    s = sqrt(s_sqr_i(:,ones(nsys+npar,1))');
    rtri = sqrt(diag(xtxi));
    rtri = rtri(:,ones(nobs+nsys,1));
    dfxhat = (x_i - x)./(s.*rtri);
    if havenans
        % Zero change in estimates if missing observations left out
        tmp = zeros(npar,length(wasnan));
        tmp(:,~wasnan) = dfxhat;
        dfxhat = tmp;
    end
    stats.(varnames{16}) = dfxhat;
end
if idx(17) % Change in fitted values. DFFIT
    dffit = h.*residuals./(1-h);
    if havenans, dffit = fixrows(dffit, wasnan); end
    stats.(varnames{17}) = dffit;
end
if idx(18) % Scaled change in fitted values. DFFITS (with studentized residuals)
    dffits = sqrt(h./(1-h)).*e_i;
    if havenans, dffits = fixrows(dffits, wasnan); end
    stats.(varnames{18}) = dffits;
end
if idx(19) %  Change in covariance. COVRATIO
    % Belsley, Kuh, and Welsch  | covratio -1 | > 3p/n -> investigate
    covr = 1 ./((((nobs-npar+ncst-1+abs(e_i).^2)./(nobs-npar+ncst)).^(npar+nsys)).*(1-h));
    if havenans, covr = fixrows(covr, wasnan); end
    stats.(varnames{19}) = covr;
end
if idx(20) %  Cook's Distance.
    d = abs(residuals).^2 .* (h./(1-h).^2)./((npar+nsys)*mse);
    if havenans, d = fixrows(d, wasnan); end
    stats.(varnames{20}) = d;
end
if idx(21) %  t Statistics.
    d = struct;
    d.xhat = xhat;
    d.se = sqrt(diag(covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    d.dfe = dfe;
    stats.(varnames{21}) = d;
end
if idx(22) %  F Statistic.
    d = struct;
    d.sse = sse;
    d.dfe = dfe;
    d.dfr = npar+nsys-1;
    d.ssr = ssr;
    d.f = (d.ssr/d.dfr)/(d.sse/d.dfe);
    d.pval = 1 - f_cdf(d.f, d.dfr, d.dfe);
    stats.(varnames{22}) = d;
end
% if idx(23) % Durbin Watson Statistic.
%     if exist('dwtest')
%         [pvalue, dw]=dwtest(residuals,Ah);
%     else
%         pvalue = 0; dw = 0;
%     end
%     d = struct;
%     d.dw = dw;
%     d.pval = pvalue;
%     stats.(varnames{23}) = d;
% end
% Oleg K. 2009 11 24 (Added HC0, HC1, HC2, HC3, HC4 and HAC statistics)
% All these computations are associated to covariance matrix estimations.
% HCi are for HCCME: Heteroskedasticity-Consistent Correlation Matrix
% Estimation (assuming correlation matrix is diagonal). HC0 is the basic
% one HC1 through HC4 are improvement w.r.t. small size samples and
% robustness to outliers.
% HAC is for Heteroskedasticity and Autocorrelation Consistent, correlation
% matrix is no more assumed to be diagonal (correlation among residuals).
% ------------------------------------------------------------------
if idx(24) % HC0
    hhat = repmat(residuals',nsys+npar,1).*Ah';
    xuux = hhat*hhat';
    d = struct;
    d.covx = xtxi*xuux*xtxi;
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{24}) = d;
end
if idx(25) % HC1
    d = struct;
    if idx(24)
        d.covx = stats.hc0.covx*(nobs+nsys)/dfe;
    else
        hhat = repmat(residuals',nsys+npar,1).*Ah';
        xuux = hhat*hhat';
        d.covx = xtxi*xuux*xtxi*(nobs+nsys)/dfe;
    end
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{25}) = d;
end
if idx(26) % HC2
    tmp = 1-h;
    xuux = zeros(nsys+npar);
    for ii = 1:(nobs+nsys)
        xuux = xuux + residuals(ii)^2/tmp(ii)*Ah(ii,:)'*Ah(ii,:);
    end
    d = struct;
    d.covx = xtxi*xuux*xtxi;
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{26}) = d;
end
if idx(27) % HC3
    tmp = (1-h).^2;
    d = struct;
    xuux = zeros(nsys+npar);
    for ii = 1:(nobs+nsys)
        xuux = xuux + residuals(ii)^2/tmp(ii)*Ah(ii,:)'*Ah(ii,:);
    end
    d.covx = xtxi*xuux*xtxi;
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{27}) = d;
end
if idx(28) % HC4
    if havenans
        if idx(10)
            h = stats.leverage;
        else
            h = fixrows(h, wasnan);
        end
    end
    tmp = (1-h).^(min(4,h/mean(h)));
    xuux = zeros(nsys+npar);
    for ii = 1:(nobs+nsys)
        xuux = xuux + residuals(ii)^2/tmp(ii)*Ah(ii,:)'*Ah(ii,:);
    end
    d = struct;
    d.covx = xtxi*xuux*xtxi;
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{28}) = d;
end
% Heteroskedasticity and Autocorrelation Consistent estimation
if idx(29) % HAC (Newey West)
    L = round(4*((nobs+nsys)/100)^(2/9));
    % L = (nobs+nsys)^.25; % as an alternative
    hhat = repmat(residuals',nsys+npar,1).*Ah';
    xuux = hhat*hhat';
    for l = 1:L
        za = hhat(:,(l+1):(nobs+nsys))*hhat(:,1:(nobs+nsys)-l)';
        w = 1 - l/(L+1);
        xuux = xuux + w*(za+za');
    end
    d = struct;
    d.covx = xtxi*xuux*xtxi;
    d.se = sqrt(diag(d.covx));
    d.t = xhat./d.se;
    d.pval = 2*(student_cdf(-abs(d.t), dfe));
    stats.(varnames{29}) = d;
end

if idx(30) % All NaNs series
    stats.empty = false;
end
if idx(31) % Rank deficient series
    stats.rankdef = false;
end
if idx(32) % Cholesky decomposition of Vinv
    stats.sqrtVinv = sqrtVinv;
end
if idx(33) % Cholesky decomposition of Winv
    stats.sqrtWinv = sqrtWinv;
end
if idx(34) % Correlation matrix from covariance matrix
    stats.corx = sqrt(diag(1./diag(covx)))*covx*sqrt(diag(1./diag(covx)));
end
% WARNING: THIERRY's COMMENT 22/6/2018
% if (idx(35)||idx(36))
%     %% Core computations with problem solving and diagnostics
%     % %% BKW Collinearity diagnostics
%     
%     [~,DD,VV] = svd(Ah./repmat(sqrt(sum(Ah.^2,1)),[size(Ah,1) 1]),'econ');
%     lambda = diag(DD(1:nsys+npar,1:nsys+npar));
%     lambda2 = lambda.*lambda;
%     VV = VV.*VV;
%     
%     phi = zeros(nsys+npar,nsys+npar);
%     for i=1:nsys+npar
%         phi(i,:) = VV(i,:)./lambda2';
%     end
%     
%     pi = zeros(nsys+npar,nsys+npar);
%     for i=1:nsys+npar
%         phik = sum(phi(i,:));
%         pi(i,:) = phi(i,:)/phik;
%     end
%     
%     lmax = lambda(1);
%     lmaxvec = lmax*ones(nsys+npar,1);
%     heta = lmaxvec./lambda;
% end
% if idx(35) % heta condition number index
%     stats.heta = heta;
% end
% if idx(36) % pi variance proportion fraction matrix
%     stats.pi = pi';
% end
end


%% Internal function for NAN removal
function vv = fixrows(v, b)
% to extend v to original length, NaNs are given by b
vv = NaN(size(b));
vv(~b) = v;
end

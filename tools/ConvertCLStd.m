function  out = ConvertCLStd(varargin)
% convert sigma into confidence level and vice versa
% for 1 and 2 dimension
p = inputParser;
p.addParameter('Mode','Sigma2CL',@(x)ismember(x,{'CL2Sigma','Sigma2CL'}));
p.addParameter('nPar',1,@(x)isfloat(x));
p.addParameter('CL',0.9,@(x)isfloat(x));
p.addParameter('Sigma',1,@(x)isfloat(x));
p.parse(varargin{:});
Mode  = p.Results.Mode;
nPar  = p.Results.nPar;
Sigma = p.Results.Sigma;
CL    = p.Results.CL;

if nPar==1
    % 1-d gaussian
    NormDist = @(x) 1./sqrt(2*pi).*exp(-0.5*x.^2);
    GetConf = @(x) arrayfun(@(a)integral(NormDist,-a,a),x);
elseif nPar==2
    % 2-d gaussian
    % this would be a square not ellipsoid/cricle
   % NormDist = @(x1,x2) 1./(2*pi).*exp(-0.5*(x1.^2+x2.^2));
   % GetConf = @(x) arrayfun(@(a)integral2(NormDist,-a,a,-a,a),x); 
    
    % polar coordinates
   % NormDist = @(r) 1./(2*pi).*exp(-0.5*r.^2);
    NormDistInt = @(r) r.*exp(-0.5*r.^2);
    GetConf = @(r) arrayfun(@(a) integral(NormDistInt,0,a),r);% arrayfun(@(a)integral(NormDistInt,0,a),r); 
  
end

switch Mode
    case 'Sigma2CL'
        out = GetConf(Sigma);
        fprintf('%.2f sigma corresponds to %.2f%% C.L. (%.0f parameter) \n',...
            Sigma,out*1e2,nPar)
    case 'CL2Sigma'
        % Get sigma as a function of confidence level
        sigmatmp = 0:0.1:20;
        [cltmp , index] = unique(GetConf(sigmatmp));
        sigmatmp = sigmatmp(index);
        GetSigma = @(cl) interp1(cltmp,sigmatmp,cl,'spline');
        
        % confidence level in as a function of sigma   
        out = GetSigma(CL);
          fprintf('%.2f%% C.L. corresponds %.2f sigma (%.0f parameter) \n',...
            CL*1e2,out,nPar)
end

end


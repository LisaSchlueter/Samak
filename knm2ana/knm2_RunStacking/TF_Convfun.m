function MaceTF_broadened = TF_Convfun(MaceTF,qU,Te,Sigma,varargin)
p = inputParser;
p.addParameter('RebinMode','Fast',@(x)ismember(x,{'Fast','Integral'}));
p.parse(varargin{:});
RebinMode = p.Results.RebinMode;

% define surplus energy
Esur = Te-qU; 
EStep = Esur(2)-Esur(1);

% add additional points at large surplus energies (does not have to be super precise, because this range is not used anyway);
EsurAdd = ((Esur(end)+EStep):EStep:(Esur(end)+EStep*200))';
Esur = [Esur;EsurAdd];
MaceTF = [MaceTF,MaceTF(end).*ones(1,numel(EsurAdd))];

% define gaussians
Gaussfun    = @(x,mu,sig) 1/(sig*sqrt(2*pi))*exp(-0.5*(x-mu).^2./sig.^2); % normalized gaussians
SingleGauss = @(x) MaceTF.*Gaussfun(x,Esur,Sigma)';                       % single gauss

% sum Gaussians and rebin
switch RebinMode
    case 'Fast'
        % just use rectangles -> fast
        % attention: doesn't work well for small Sigma < 0.1 eV
        SumGauss = @(x) sum(SingleGauss(x),2).*EStep;
        MaceTF_broadened = SumGauss(Esur'); % superposition, this isn't a proper integration
    case 'Integral'
        % proper integration -> slow
        SumGauss = @(x) sum(SingleGauss(x),2);
        MaceTF_broadened = arrayfun(@(a,b) integral(SumGauss,a,b,'ArrayValued',1) ,Esur'-EStep/2,Esur'+EStep/2);
end

MaceTF_broadened = MaceTF_broadened(1:numel(Te));   % cut energy range again to initial range
end
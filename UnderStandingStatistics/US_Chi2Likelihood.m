Mode = 'Gauss2';%'Gauss';
mu    = 0;
sigma = 1;

switch Mode
    case 'Gauss'
        % draw from distribution
        xRand = mu+sigma.*randn(1e4,1); % -> fit values
         Model1 = @(x)(mu-x).^2./sigma^2;
         Model2 = @(x)(mu-x).^2./sigma^2;
         Model = @(x)(mu-x).^2./sigma^2;
         Y = Model(mu);
    case 'Gauss2'
        nRand  = 3e4;
        sigma2 = sigma.*1.2;
        xRand1 = mu+sigma.*randn(nRand,1); % -> fit values
        xRand2 = mu+sigma2.*randn(nRand,1); % -> fit values
        xRand  = [xRand1(xRand1>=mu);xRand2(xRand2<mu)];
        
        Model1 = @(x) (mu-x).^2./sigma^2;
        Model2 = @(x) (mu-x).^2./(sigma2)^2;
        Model  = @(x) (x>=mu).*Model1(x)+(x<mu).*Model2(x);%@(x) (x>=mu).*(mu-x).^2./sigma^2+(x<mu).*(mu-x).^2./(sigma2)^2;
        Y      = Model(mu);
end
close all;
x = linspace(min(xRand),max(xRand),1e3);
chi2 = chi2func(Model(x),Model(mu),sqrt(Model(x)));
%%
%Probtmp = @(y)  exp(-0.5*chi2func(Model(y),Model(mu),sqrt(Model(y))))./simpsons(x,exp(-0.5.*chi2));
Norm1 = 2*simpsons(x(x>=mu),exp(-0.5.*chi2(x>=mu)));
Norm2 = 2*simpsons(x(x<mu),exp(-0.5.*chi2(x<mu)));
Probtmp=@(y) (y>=mu).* exp(-0.5*chi2func(Model1(y),Model1(mu),sqrt(Model1(y))))./Norm1+...
              (y<mu).* exp(-0.5*chi2func(Model2(y),Model2(mu),sqrt(Model2(y))))./Norm2;
          
%Probtmp = @(y)  exp(-0.5*chi2func(Model(y),Model(mu),sqrt(Model(y))))./simpsons(x,exp(-0.5.*chi2));
simpsons(x(x<=mu),Probtmp(x(x<=mu)))
simpsons(x(x>mu),Probtmp(x(x>mu)))
Prob = Probtmp(x);
%%
%Prob = exp(-0.5.*chi2)./simpsons(x,exp(-0.5.*chi2));
% 
% plot(x,chi2)
% grid on;
GetFigure
h1 = histogram(xRand,'Normalization','pdf');
hold on;
p1 = plot(x,Prob,'LineWidth',2);

function chi2 = chi2func(x,mu,sigmaModel)
chi2 = (x-mu).^2./sigmaModel.^2;
end


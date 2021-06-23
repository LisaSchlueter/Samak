

muTrue    = 0.2;
sigmaTrue = 1;

%% 1. find [x1,x2] with MC
nToyMC = 5e4;
xMC = sort(muTrue+sigmaTrue.*randn(nToyMC,1));
chi2TrueMC = ComputeChi2(xMC,muTrue,sigmaTrue);
chi2BestMC = ComputeChi2(xMC,[zeros(numel(xMC(xMC<0)),1);xMC(xMC>=0)],sigmaTrue);
Deltachi2MC = chi2TrueMC-chi2BestMC;

% GetFigure;
% plot(xMC,Deltachi2MC)

Chi2Crit_v = 0.8:0.01:2.5;
CumProb = sum(Deltachi2MC<=Chi2Crit_v)./nToyMC;
DeltaChi2CritMC = interp1(CumProb,Chi2Crit_v,0.9,'spline');
x1_MC = interp1(Deltachi2MC(xMC<0),xMC(xMC<0),DeltaChi2CritMC,'spline');
x2_MC = interp1(Deltachi2MC(xMC>0),xMC(xMC>0),DeltaChi2CritMC,'spline');

%% 2. find [x1,x2] with Probability trick
nToy = 300;
xToy = linspace(min(xMC),max(xMC),nToy)';
chi2True = ComputeChi2(xToy,muTrue,sigmaTrue);
chi2Best = ComputeChi2(xToy,[zeros(numel(xToy(xToy<0)),1);xToy(xToy>=0)],sigmaTrue);
Deltachi2 = chi2True-chi2Best;

Prob_tmp = exp(-0.5*chi2True);
Prob = Prob_tmp./simpsons(xToy,Prob_tmp);
CumProb = GetCDF(xToy,Prob);

% GetFigure;
% plot(xToy,Prob); hold on;
% plot(xToy,CumProb)

interpStyle = 'spline';

a =  @(x) (interp1(xToy,CumProb,x(2),interpStyle)-interp1(xToy,CumProb,x(1),interpStyle))-0.9;
b =  @(x)  (interp1(xToy,Deltachi2,x(1),interpStyle)-interp1(xToy,Deltachi2,x(2),interpStyle));
fun = @(x) [a(x),b(x)];
x0 = [-1,1];
options = optimoptions(@fsolve,'Display','off');
x1x2 = fsolve(fun,x0,options);
x1_Toy = x1x2(1);
x2_Toy = x1x2(2);
DeltaChi2CritToy = interp1(xToy,Deltachi2,x1_Toy,interpStyle);

%% results
fprintf('MC : x1 = %.3f , x2 = %.3f , chi2crit = %.3f \n',x1_MC,x2_MC,DeltaChi2CritMC)
fprintf('Toy: x1 = %.3f , x2 = %.3f , chi2crit = %.3f\n',x1_Toy,x2_Toy,DeltaChi2CritToy)

                    %%
function chi2 = ComputeChi2(x,mu,sigma)

chi2 = (x-mu).^2./sigma.^2;

end
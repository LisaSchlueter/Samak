% Dummy script to understand why KNM-2 FSD uncertainty is different from Saenz 2000
sigma = 0.5;
mu = 2.5;
DistFun = @(x) 1./sqrt(2*pi*sigma^2).*exp(-0.5*(x-mu).^2./sigma^2);

Emean1 = 0.1:0.2:4.9;
Prob1 = arrayfun(@(a,b) integral(DistFun,a,b,'AbsTol',1e-9),Emean1-0.1,Emean1+0.1);

Emean2 = 0.05:0.1:4.95;
Prob2 = arrayfun(@(a,b) integral(DistFun,a,b,'AbsTol',1e-9),Emean2-0.05,Emean2+0.05);

Emean3 = 0.025:0.05:4.75;
Prob3 = arrayfun(@(a,b) integral(DistFun,a,b,'AbsTol',1e-9),Emean3-0.025,Emean3+0.025);

Var1 = var(Emean1,Prob1);
Var2 = var(Emean2,Prob2);
Var3 = var(Emean3,Prob3);

nSamples = 1e3;
% randomize
Prob2_rand = Prob2.*(1+randn(nSamples,numel(Prob2)).*0.05);
Prob1_rand = Prob1.*(1+randn(nSamples,numel(Prob1)).*0.05);
Prob3_rand = Prob3.*(1+randn(nSamples,numel(Prob3)).*0.05);

close all
boundedline(Emean1,mean(Prob1_rand),std(Prob1_rand));
hold on;
boundedline(Emean2,mean(Prob2_rand),std(Prob2_rand));
boundedline(Emean3,mean(Prob3_rand),std(Prob3_rand));

Var1_samples = zeros(nSamples,1);
Var2_samples = zeros(nSamples,1);
Var3_samples = zeros(nSamples,1);
for i=1:nSamples
    Var1_samples(i) = var(Emean1,Prob1_rand(i,:));  
    Var2_samples(i) = var(Emean2,Prob2_rand(i,:)');  
    Var3_samples(i) = var(Emean3,Prob3_rand(i,:)');  
end


fprintf('rel err variance %.2g binning = %.2f%% \n',Emean1(2)-Emean1(1),std(Var1_samples)./Var1*100);
fprintf('rel err variance %.2g binning = %.2f%% \n',Emean2(2)-Emean2(1),std(Var2_samples)./Var2*100);
fprintf('rel err variance %.2g binning = %.2f%% \n',Emean3(2)-Emean3(1),std(Var3_samples)./Var3*100);

%RelErrVar2 = std(Var2_samples)./Var2*100;

return
% Dummy test
E    = [0.5,1.5,2.5];
Prob = [0.1,0.2,0.1];
Prob = Prob./simpsons(E,Prob);
bar(E,Prob)


E2    = 0.25:0.5:2.75;
Prob2 = [0.05,0.05,0.1,0.1,0.05,0.05];
Prob2 = Prob2./simpsons(E2,Prob2);
hold on;
bar(E2,Prob2);

var(E,Prob);
var(E2,Prob2);

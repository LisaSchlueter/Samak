% Toy MC to understand difference marginal vs profile likelihood

% 2-dimensional scenario
x_true = 2;
nSamples = 5e3;
r = sort(poissrnd(x_true,nSamples,1));

x = -5:0.1:x_true+5*x_true;
chi2 = (x-x_true).^2./x_true;
prob = exp(-0.5.*chi2)./simpsons(x,exp(-0.5.*chi2));


%% plot
close all;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.35]);
s1 = subplot(1,2,1);
h1 = histogram(r,'Normalization','pdf','FaceColor',rgb('DodgerBlue'));
hold on;
p1 =plot(x,prob,'LineWidth',2);
xlabel('x');
ylabel('pdf');
PrettyFigureFormat;
xlim(x_true+[-3*sqrt(x_true),3*sqrt(x_true)]);

s2 = subplot(1,2,2);
plot(linspace(min(x),max(x),10),0.5.*ones(10,1),'k:','LineWidth',2);
hold on;
plot(x_true.*ones(10,1),linspace(-0.05,1.05,10),'k:','LineWidth',2);
plot(r,cumsum(r)./max(cumsum(r)),'LineWidth',2,'Color',h1.FaceColor);
plot(x,GetCDF(x,prob),'LineWidth',2,'Color',p1.Color);
xlabel('x');
ylabel('cdf');
PrettyFigureFormat;
ylim([-0.05 1.05])
xlim(x_true+[-3*sqrt(x_true),3*sqrt(x_true)]);
linkaxes([s1,s2],'x');

% Read Fabian Data
F=importdata('28032019_ProKatrin_E0ShiftsRingwise_withoutErrors.txt');
R=1:10;
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 1000])
for i=1:1:34
hold on
scatter(R,F(R,i+1),'d','filled')
end
xlabel('ring');
ylabel('\Delta U');
hold off
PrettyFigureFormat
set(gca,'FontSize',24);


%% Interpolation
G=F(R,2:end)+E090_S(R)';
ring = [R];
bias = [-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 1.1 1.2]; 
[x,y] = meshgrid(ring,bias);
eq = @(xq,yq) griddata(ring,bias,G',xq,yq);

figure(777)
set(gcf, 'Position',  [100, 100, 1000, 1000])
for i=bias
hold on
plot(ring,eq(ring,i))
end
hold on
h=plot(ring,eq(ring,0.175),'Color','Red','LineWidth',5);
d=errorbar(ring,E090_S(R),E090Err_S(R),'ks','MarkerSize',10,'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',2) ;
hold off
xlabel('ring');
ylabel('E_0 (eV)');
hold off
legend([h,d],'Current RW Settings: +175 mV','Samak Data');legend boxoff;
PrettyFigureFormat
set(gca,'FontSize',24);


%% Samak
E090_S    = [18573.181        18573.27       18573.399       18573.553       18573.444       18573.479       18573.468       18573.571       18573.561       18573.704       18573.547       18573.677];
E090Err_S = [0.134       0.077       0.077       0.076       0.077       0.077       0.077       0.077       0.085       0.089       0.096       0.154];
E040_S    = [18573.245       18573.317       18573.516       18573.631       18573.558       18573.557       18573.476       18573.797       18573.617       18573.697       18573.758        18573.59];
E040Err_S = [0.237       0.137       0.138       0.136       0.138       0.137       0.137       0.138       0.153       0.161       0.172       0.281];
E090m_S     = wmean(E090_S,1./E090Err_S.^2);
E040m_S     = wmean(E040_S,1./E040Err_S.^2);

%% KaFit DeepScans
WQ=importdata('Wongkook.txt')
E0K_m1V    = WQ(:,2);
E0ErrK_m1V = WQ(:,3);
E0K_m02V    = WQ(:,4);
E0ErrK_m02V = WQ(:,5);
E0K_00V    = WQ(:,6);
E0ErrK_00V = WQ(:,7);
E0K_p03V    = WQ(:,8);
E0ErrK_p03V = WQ(:,9);

%% check Plot
myring = [R];

figure(100)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,0.)-1,'LineWidth',5);
hold on
d=errorbar(myring,E0K_00V(myring),E0ErrK_00V(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model','Wongkook Fit'); legend boxoff;
title('Deep E0 scans- RW Bias: 0V')
PrettyFigureFormat

figure(101)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,-1)-1,'LineWidth',5);
hold on
d=errorbar(myring,E0K_m1V(myring),E0ErrK_m1V(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model','Wongkook Fit'); legend boxoff;
title('Deep E0 scans- RW Bias: -1V')
PrettyFigureFormat

figure(102)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,-0.2)-1,'LineWidth',5);
hold on
d=errorbar(myring,E0K_m02V(myring),E0ErrK_m02V(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model','Wongkook Fit'); legend boxoff;
title('Deep E0 scans- RW Bias: -0.2V')
PrettyFigureFormat

figure(103)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,0.3)-1,'LineWidth',5);
hold on
d=errorbar(myring,E0K_p03V(myring),E0ErrK_p03V(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
legend([m,d],'Ferenc Model','Wongkook Fit'); legend boxoff;
hold off
title('Deep E0 scans- RW Bias: 0.3V')
PrettyFigureFormat

%%
%% Wongkook 300mV
WQ=importdata('Wongkooq300mV.txt')
E0K_300mV     = WQ(:,2);
E0ErrK_300mV  = WQ(:,3);

%% check Plot
myring = [R];
figure(200)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,0.),'LineWidth',5);
hold on
d=errorbar(myring,E0K_300mV(myring)+0.35,E0ErrK_300mV(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('DarkGreen'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model','Wongkook Fit','Location','SouthEast'); legend boxoff;
title('Regular Beta Scans scans- RW Bias: 300 mV')
PrettyFigureFormat

%%
figure(2)
set(gcf, 'Position',  [100, 100, 1200, 800])
for i=1:1:34
hold on
h(i)=plot(R,F(R,i+1)+E090_S(R)'-E090m_S,'LineWidth',5)
end
xlabel('ring');
ylabel('Corrected E_0 - 18573.48 (eV)');
%legend([h(1) h(2) h(3) h(4) h(5)],'-1V','-0.2V','0V','+0.3V','+1V')
hold off
title('Predicted Ring-Dependent Endpoint');
PrettyFigureFormat
set(gca,'FontSize',20);


%% Estimator 1 : standard Deviation over rings
myring = [R];
Counter=0;
for i = bias
    Counter=Counter+1;
    estimator1(Counter) = (std(eq(myring,i)));
end
figure(3)
set(gcf, 'Position',  [100, 100, 500, 800])
plot(mybias,estimator1,'LineWidth',5)
ylabel('Minimum STD Estimator');
xlabel('Rear Wall bias (V)');
PrettyFigureFormat;

%% Estimator 2 : fit slope
myring = [R];
Counter=0;
for i = bias
    Counter=Counter+1;
     P(Counter,:) = polyfit(myring,eq(myring,i),1);
end
estimator2 = abs(P(:,1));
figure(4)
set(gcf, 'Position',  [100, 100, 500, 800])
plot(mybias,estimator2,'LineWidth',5)
ylabel('Minimum Slope Estimator');
xlabel('Rear Wall bias (V)');
PrettyFigureFormat;

%% check Plot xxx
myring = [R];
figure(5)
set(gcf, 'Position',  [100, 100, 1000, 600])
d=plot(myring,eq(myring,0.175),'LineWidth',5)
hold on
errorbar(myring,E090_S(R),E090Err_S(R),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',5) ;
p1=plot(myring,eq(myring,-0.175)-0.275,'LineWidth',8,'Color',rgb('DarkGreen'));
p2=plot(myring,eq(myring,-0.330)-0.4,'LineWidth',8,'Color',rgb('Orange'));
ylabel('Endpoint (eV)');
xlabel('Ring');
hold off
legend([d,p1,p2],'Actual RW 175 mV','Prediction -175 mV','Prediction -350 mV','Location','SouthEast');  legend boxoff;
PrettyFigureFormat;
grid on
set(gca,'FontSize',24);



%% Summary Plot
figure(6)
set(gcf, 'Position',  [100, 100, 1200, 800])
xlabel('Rear Wall bias (V)');
yyaxis left
plot(mybias,smooth(estimator2,'rlowess'),'LineWidth',5)
ylabel('Minimum Slope Estimator');
yyaxis right
plot(mybias,smooth(estimator1,'rlowess'),'LineWidth',5)
ylabel('Minimium STD Estimator');
PrettyFigureFormat
title('RW Bias Voltage Optimization: 2:11 Rings');
set(gca,'FontSize',24);

%% Summary Plot
figure(7)
set(gcf, 'Position',  [100, 100, 1000, 800])
plot(mybias,smooth(estimator2.*estimator1','rlowess'),'LineWidth',15,'Color',rgb('SteelBlue'));
%plot(mybias,(estimator2.*estimator1'),'LineWidth',5,'Color',rgb('SteelBlue'));
xlabel('Rear Wall bias (V)');
ylabel('Minimum Slope x Minimum STD Estimator');
PrettyFigureFormat
set(gca,'FontSize',24);
xlim([-0.4 -0.1]);
grid on

%%
figure(15)
ring = [R];
bias = [-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 1.1 1.2];
%[x,y] = meshgrid(ring,bias);
ci=0;cj=0;
for i=ring
    ci=ci+1;
    cj=0;
    for j=bias
        cj=cj+1;
        mye0(ci,cj)=eq(i,j);
    end
end
ribbon(mye0);
xlabel('Rear Wall bias (V)');
ylabel('Ring');
zlabel('Predicted Endpoint (eV)');
PrettyFigureFormat



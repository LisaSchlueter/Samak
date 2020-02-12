%% Reference Values KaFit/Samak @-150mV
E090_Wm150    = [ 18573.66596809, 18573.71700657, 18573.68761071, 18573.6398121, 18573.67854766, 18573.76886288, 18573.69531948, 18573.70395726, 18573.70512517, 18573.6945384, 18573.79023627, 18573.74122678 ];
E090Err_Ref150= [ 0.06115680596795, 0.03515980143713, 0.03517945954824, 0.03511380656838, 0.03515706216064, 0.03521480558135, 0.03520002116202, 0.03527094100592, 0.04067956130441, 0.04061295230847, 0.04326569072667, 0.07089417824411];
E090_Ref150   = E090_Wm150;

%% Fit Values Samak @-50mV
E090_Sm50     = 18573.61 + [-0.0570    0.0801   -0.0781   -0.0733   -0.0269   -0.0927   -0.0587   -0.0722   -0.0794   -0.0379   -0.0940   -0.0099];


%% Read Fabian Data
%  gives the relative endpoint shift per ring and per rear wall voltage. 
%  The data is RELATIVE TO U_RW=-0.150V  
F=importdata('19092019_ProKatrin_ExcludeLastScan_E0ShiftsRingwise_withoutErrors.txt');
R=1:12;
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 1000])
%for i=1:1:40
select=5:1:30;
%select=1:1:40;
for i=select
hold on
scatter(R,F(R,i+1),'d','filled')
end
colorbar
xlabel('ring');
ylabel('\Delta U (V) - Change of Reference Potential');
hold off
PrettyFigureFormat
set(gca,'FontSize',24);


%% Interpolation
%  eq(RING,BIAS) --> Endpoint Shift
%G=F(R,2:end)+E090_Wm150(R)';
G=+F(R,2+4:2+29)+E090_Ref150(R)';
ring = [R];
bias = [-3.0 -2.0 -1.5 -1.2 -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2.0 3.0 4.0]; 
bias = bias(select);
[x,y] = meshgrid(ring,bias);
eq = @(xq,yq) griddata(ring,bias,G',xq,yq);
figure(777)
set(gcf, 'Position',  [100, 100, 1000, 1000])
for i=bias
hold on
plot(ring,eq(ring,i))
end
hold on
h=plot(ring,eq(ring,-0.150),'Color','Red','LineWidth',5);
d=errorbar(ring,E090_Ref150(R),E090Err_Ref150(R),'ks','MarkerSize',10,'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',2) ;
hold off
xlabel('ring');
ylabel('E_0 (eV)');
hold off
legend([h,d],'Reference Masurement at RW = -150 mV','Samak Data','Location','SouthEast');legend boxoff;
PrettyFigureFormat
set(gca,'FontSize',24);

%% Estimator 1 : standard Deviation over rings
myring = [R];
Counter=0;
for i = bias
    Counter=Counter+1;
    estimator1(Counter) = (std(eq(myring,i)));
end
figure(3)
set(gcf, 'Position',  [100, 100, 500, 800])
plot(bias,estimator1,'LineWidth',5)
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
plot(bias,estimator2,'LineWidth',5)
ylabel('Minimum Slope Estimator');
xlabel('Rear Wall bias (V)');
PrettyFigureFormat;


%% Summary Plot
figure(6)
set(gcf, 'Position',  [100, 100, 1200, 800])
xlabel('Rear Wall bias (V)');
yyaxis left
plot(bias,smooth(estimator2,'rlowess'),'LineWidth',5)
ylabel('Minimum Slope Estimator');
yyaxis right
plot(bias,smooth(estimator1,'rlowess'),'LineWidth',5)
ylabel('Minimium STD Estimator');
PrettyFigureFormat
title('RW Bias Voltage Optimization: 2:11 Rings');
set(gca,'FontSize',24);

%% Summary Plot
figure(7)
set(gcf, 'Position',  [100, 100, 1000, 800])
plot(bias,smooth(estimator2.*estimator1','rlowess'),'LineWidth',15,'Color',rgb('SteelBlue'));
%plot(mybias,(estimator2.*estimator1'),'LineWidth',5,'Color',rgb('SteelBlue'));
xlabel('Rear Wall bias (V)');
ylabel('Minimum Slope x Minimum STD Estimator');
PrettyFigureFormat
set(gca,'FontSize',24);
xlim([-0.5 0.3]);
grid on

%%
figure(15)
ring = [R];
bias = [-3.0 -2.0 -1.5 -1.2 -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2.0 3.0 4.0]; 
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

%%
figure(101)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,-0.05),'LineWidth',5);
hold on
d=errorbar(myring,E090_Sm50(myring),E090Err_Ref150(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model RW@-50mV','Samak Fit RW@-50mV'); legend boxoff;
title('E_0 Ring-wise -  RW@-50mV')
PrettyFigureFormat


%%
figure(102)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,-0.150),'LineWidth',5);
hold on
d=errorbar(myring,E090_Wm150(myring),E090Err_Ref150(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model RW@-150mV','KaFit Fit RW@-150mV'); legend boxoff;
title('E_0 Ring-wise -  RW@-150mV')
PrettyFigureFormat

%%
figure(102)
set(gcf, 'Position',  [100, 100, 1000, 300])
m=plot(myring,eq(myring,-0.05),'LineWidth',5);
hold on
d=errorbar(myring,E090_Sm50(myring),E090Err_Ref150(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
m2=plot(myring,eq(myring,-0.2),'LineWidth',5);
hold off
legend([m,d,m2],'Ferenc Model RW@-50mV','Samak Fit RW@-50mV','Optimum RW@-200mV'); legend boxoff;
PrettyFigureFormat
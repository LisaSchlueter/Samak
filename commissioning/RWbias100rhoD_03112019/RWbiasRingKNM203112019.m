%% Reference Values KaFit/Samak @-8mV
E090_Wm150    = 18575.7 + [ 0.0534    0.0496    0.0225    0.0514    0.0242   -0.0135   -0.0421   -0.0153   -0.0238   -0.0009   -0.0699    0.0483 ];
E090Err_Ref150= [ 0.0486    0.0282    0.0283    0.0284    0.0287    0.0289    0.0289    0.0291    0.0321    0.0340    0.0390    0.0591];
E090_Ref150   = E090_Wm150;

%% Fit Values Samak @-50mV
E090_Sm50     = E090_Wm150;

%% Read Fabian Data
%  gives the relative endpoint shift per ring and per rear wall voltage. 
%  The data is RELATIVE TO U_RW=-0.150V  
%F=importdata('01112019_ProKatrin_UMS1_E0ShiftsRingwise_withoutErrors.txt');
F=importdata('03112019_ProKatrin_UMS2_E0ShiftsRingwise_withoutErrors.txt');
R=1:11;
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 1000])
select=1:1:16;
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
G=F(R,2:end)+E090_Wm150(R)';
ring = [R];
bias = [-0.2 -0.15 -0.1 -0.05 -0.008 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6];
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
h=plot(ring,eq(ring,-0.008),'Color','Red','LineWidth',5);
d=errorbar(ring,E090_Ref150(R),E090Err_Ref150(R),'ks','MarkerSize',10,'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',2) ;
hold off
xlabel('ring');
ylabel('E_0 (eV)');
hold off
legend([h,d],'Reference Masurement at RW = -8 mV','Samak Data','Location','SouthEast');legend boxoff;
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
bias = [-0.2 -0.15 -0.1 -0.05 -0.008 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6];
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
m=plot(myring,eq(myring,-0.008),'LineWidth',5);
hold on
d=errorbar(myring,E090_Sm50(myring),E090Err_Ref150(myring),'ks','MarkerSize',20,'MarkerFaceColor',rgb('Amethyst'),'LineWidth',2) ;
hold off
legend([m,d],'Ferenc Model RW@-8mV','Samak Fit RW@-8mV'); legend boxoff;
title('E_0 Ring-wise -  RW@-8mV')
PrettyFigureFormat


filename = 'RAA2011_95a.txt';
d = importdata(filename);
filename2 = 'RAA2011_95b.txt';
d2 = importdata(filename2);
filename3 = 'RAA2011_95cnew.txt';
d3 = importdata(filename3);

close all;
GetFigure;
plot(d(:,1),d(:,2),'-','LineWidth',3);
hold on
plot(d(1:15:end,1),smooth(d(1:15:end,2),5),'-.','LineWidth',2,'MarkerSize',10);

plot(d2(:,1),d2(:,2),'-','LineWidth',3);
hold on
plot(d2(1:15:end,1),smooth(d2(1:15:end,2),5),'-.','LineWidth',2,'MarkerSize',10);

plot(d3(:,1),d3(:,2),'-','LineWidth',3);
hold on
plot(d3(1:15:end,1),smooth(d3(1:15:end,2),5),'-.','LineWidth',2,'MarkerSize',10);
set(gca,'XScale','log');
set(gca,'YScale','log');
ylim([0.1 3])

path = [getenv('SamakPath'),'inputs/FSD/'];

Iso = 'T2';
f0  = [path,'FSD_KNM2_',Iso,'.txt'];
f1  = [path,'FSD_KNM2_',Iso,'_0p1eV.txt'];
f2  = [path,'FSD_KNM2_',Iso,'_0p5eV.txt'];

d0 = importdata(f0);
d1 = importdata(f1);
d2 = importdata(f2);

GetFigure;

p0 = plot(d0(:,1),d0(:,2),'-','LineWidth',2);
hold on;
p1 = plot(d1(:,1),d1(:,2),'-.','LineWidth',2);
p2 = plot(d2(:,1),d2(:,2),':','LineWidth',2);
PrettyFigureFormat;
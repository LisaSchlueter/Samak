 
% plot and compare some KNM2 FSD
Iso = 'T2';

% load FSDs from file
path = [getenv('SamakPath'),'inputs/FSD/'];
d           = importdata([path,'/FSD_KNM2_',Iso,'.txt']);
d_0p1       = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV.txt']);
d_0p1_cut40 = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV_cut40eV.txt']);
d_0p1_cut50 = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV_cut50eV.txt']);

%% plot

GetFigure;
p1 = plot(d(:,1),d(:,2),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(d_0p1(:,1),d_0p1(:,2),'-.','LineWidth',2,'Color',rgb('Orange'));
p3 = plot(d_0p1_cut40(:,1),d_0p1_cut40(:,2),':','LineWidth',2,'Color',rgb('ForestGreen'));
hold on;
p4 = plot(d_0p1_cut50(:,1),d_0p1_cut50(:,2),'--','LineWidth',2,'Color',rgb('FireBrick'));

PrettyFigureFormat;
leg = legend([p1,p2,p3,p4],...
    sprintf('KNM2 default (%.0f bins)',numel(d(:,1))),...
    sprintf('KNM2 0.1eV (%.0f bins)',numel(d_0p1(:,1))),...
    sprintf('KNM2 0.1eV - cut 40 eV (%.0f bins)',numel(d_0p1_cut40(:,1))),...
    sprintf('KNM2 0.1eV - cut 50 eV (%.0f bins)',numel(d_0p1_cut50(:,1))),...
    'Location','southeast');
PrettyLegendFormat(leg);
xlabel('Energy (eV)');
xlim([-3 60]);
set(gca,'YScale','log');
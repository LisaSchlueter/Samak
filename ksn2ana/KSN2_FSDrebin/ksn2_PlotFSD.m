 
% plot and compare some KNM2 FSD
Iso = 'T2';

% load FSDs from file
path = [getenv('SamakPath'),'inputs/FSD/'];
d           = importdata([path,'/FSD_KNM2_',Iso,'.txt']);
d_0p1       = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV.txt']);
d_0p5       = importdata([path,'/FSD_KNM2_',Iso,'_0p5eV.txt']);
%d_0p1_cut40 = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV_cut40eV.txt']);
% d_0p1_cut50 = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV_cut50eV.txt']);

%% plot

GetFigure;
p1 = plot(d(:,1),d(:,2),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(d_0p1(:,1),d_0p1(:,2),'-.','LineWidth',2,'Color',rgb('Orange'));
p3 = plot(d_0p5(:,1),d_0p5(:,2),':','LineWidth',2,'Color',rgb('LimeGreen'));
% hold on;
% p4 = plot(d_0p1_cut50(:,1),d_0p1_cut50(:,2),'--','LineWidth',2,'Color',rgb('FireBrick'));

PrettyFigureFormat;
leg = legend([p1,p2,p3],...
    sprintf('KNM2 default (%.0f bins)',numel(d(:,1))),...
    sprintf('KNM2 0.1eV (%.0f bins)',numel(d_0p1(:,1))),...
    sprintf('KNM2 0.5eV (%.0f bins)',numel(d_0p5(:,1))),...
    'Location','southeast');
%    sprintf('KNM2 0.1eV - cut 40 eV (%.0f bins)',numel(d_0p1_cut40(:,1))),...
%    sprintf('KNM2 0.1eV - cut 50 eV (%.0f bins)',numel(d_0p1_cut50(:,1))),...
PrettyLegendFormat(leg);
xlabel('Energy (eV)');
xlim([-3 60]);
set(gca,'YScale','log');

%% plot binning
GetFigure;

stairs(d(2:end,1),diff(d(:,1)),'LineWidth',2)
hold on;
stairs(d_0p1(2:end,1),diff(d_0p1(:,1)),'LineWidth',2)
grid on;
xlim([-0.5 55]);
ylim([0,1.1]);
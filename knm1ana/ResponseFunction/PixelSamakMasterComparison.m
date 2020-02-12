% Comparison of Master & Samak RF 
% Pixel-Wise
close all

load('SamakRF.mat');
load('MasterRF.mat');

%% Definition
E0    = 18574;
qu    = [18373  18483  18488  18493  18498  18503  18508  18513  18518  18523  18528  18533  18535  18537  18539  18541  18543  18545  18547  18549  18551  18553  18555  18557  18559  18561  18562  18563  18564  18565  18566  18567  18569  18571  18573  18578  18583  18593  18603  18623];
qUbin = 12;
pixel = 90;

%% SuperImpose Master / Samak
fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
r=stairs(MasterTe(:,pixel)-E0,MasterRF(qUbin,:,pixel),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
hold on
s=stairs(SamakTe-E0,SamakRF(qUbin,:,pixel),'-','Color',rgb('IndianRed'),'LineWidth',2);
hold off
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',18);
ylabel(sprintf('transmission pixel%.0f',pixel),'FontSize',18);
legend([r s],'Master RF','Samak RF','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',18);

% SuperImpose Master / Samak - Interpolated
fMasterRF = @(te) interp1(MasterTe(:,pixel),MasterRF(qUbin,:,pixel),te);
fSamakRF = @(te) interp1(SamakTe,SamakRF(qUbin,:,pixel),te);
te=18320:0.001:18600

fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
r=stairs(te-E0,fMasterRF(te),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
hold on
s=stairs(te-E0,fSamakRF(te),'-','Color',rgb('IndianRed'),'LineWidth',2);
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',18);
ylabel(sprintf('transmission pixel%.0f',pixel),'FontSize',18);
legend([r s],'Master RF - Interpolated','Samak RF - Interpolated','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',18);

fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
r=stairs(te,fSamakRF(te)-fMasterRF(te),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',18);
ylabel(sprintf('transmission difference pixel%.0f',pixel),'FontSize',18);
legend([r],'Samak - Master RF','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',18);
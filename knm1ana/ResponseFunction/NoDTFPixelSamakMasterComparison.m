% Comparison of Master & Samak RF 
% NO DETAILED TRANSMISSION
% Pixel-Wise

% Samak
a=importdata('SamakRFqU18545.mat');
SamakRF(:,:,1) = a.SamakRF(:,:);
SamakTe        = a.SamakTe;

% Fitrium
%fitrium=importdata('response_fitrium_p0_18545.txt');

%% MasterNodtfRF
load('MasterNodtfRF.mat');
MasterTe                     = MasterNodtfTe;
for i=1:5
MasterRF(i,:,:)              = MasterNodtfRF;
end


%% Definition
E0    = 18574;
qu    = [18373 18523 18545 18573 18623];
qUbin = 4;
pixel = 1;

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
%f=stairs(fitrium(:,1)-E0+18600 + -18598.0384,fitrium(:,2),'-','Color',rgb('DarkGreen'),'LineWidth',2);
hold off
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',24);
ylabel(sprintf('transmission pixel%.0f',pixel),'FontSize',24);
xlim([qu(qUbin)-5-E0 SamakTe(end)-E0]);
%legend([r s f],'Master RF','Samak RF','Fitrium RF','Location','NorthWest');
legend([r s],'Master RF','Samak RF','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',24);

%% SuperImpose Master / Samak - Interpolated
fMasterRF = @(te) interp1(MasterTe(:,pixel),MasterRF(qUbin,:,pixel),te);
fSamakRF = @(te) interp1(SamakTe,SamakRF(qUbin,:,pixel),te);
%fFitriumRF = @(te) interp1(fitrium(:,1)+18600 -18598.0384,fitrium(:,2),te);
te=qu(qUbin)-5:0.1:18600;

fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';

s1=subplot(2,1,1)
r=stairs(te-E0,fMasterRF(te),'-','Color',rgb('DarkSlateGray'),'LineWidth',4);
hold on
s=stairs(te-E0,fSamakRF(te),'-','Color',rgb('IndianRed'),'LineWidth',2);
%f=plot(te-E0,fFitriumRF(te),'-','Color',rgb('DarkGreen'),'LineWidth',2);
hold off
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',20);
ylabel(sprintf('transmission pixel%.0f',pixel),'FontSize',18);
%legend([r s f],'Master RF - Interpolated','Samak RF - Interpolated','Fitrium RF','Location','NorthWest');
legend([r s],'Master RF','Samak RF','Location','Best'); legend('boxoff');
grid on
xlim([qu(qUbin)-E0 qu(qUbin)+3-E0]);
PrettyFigureFormat
set(gca,'FontSize',18);

s2=subplot(2,1,2)

a.FontSize=20;a.FontWeight='bold';
r1=stairs(te-E0,(fSamakRF(te)-fMasterRF(te)),'-','Color',rgb('DarkSlateGray'),'LineWidth',4);
hold on
%r2=stairs(te-E0,(-fSamakRF(te)+fFitriumRF(te)),'-','Color',rgb('IndianRed'),'LineWidth',2);
%r3=stairs(te-E0,(+fMasterRF(te)-fFitriumRF(te)),'-','Color',rgb('CadetBlue'),'LineWidth',2);
hold off
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e -%.1f (eV)',E0),'FontSize',20);
ylabel(sprintf('transmission difference pixel%.0f',pixel),'FontSize',18);
%legend([r1 r2 r3],'Samak - Master RF','Samak - Fitrium RF','Fitrium - Master RF','Location','NorthWest');
legend([r1 ],'Samak - Master RF','Location','Best');legend('boxoff');
grid on
xlim([qu(qUbin)-E0 qu(qUbin)+3-E0]);
PrettyFigureFormat
set(gca,'FontSize',18);
legend([r1 ],'Samak - Master RF','Location','Best','FontSize',24);legend('boxoff');

linkaxes([s1,s2],'x');

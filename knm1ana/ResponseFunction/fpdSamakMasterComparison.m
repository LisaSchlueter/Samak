% Comparison of Master & Samak RF 
% Averaged over a given Pixel List
ComputeSamak  = 'ON';
GetMaster     = 'ON';

if strcmp(GetMaster,'ON') 
% Get Master RFs
[MasterTe MasterRFfpd] = fpdAverageMasterRF();
end
if strcmp(ComputeSamak,'ON') 
% Get Samak RFs
!rm ../..//inputs/ResponseFunction/samakRF_Uniform_1.1e+17cm2_NIS7_Bm4.23T_Bs2.52T_Ba6.3036G_Temin18373.000_Temax18623.000_Bin100_KatrinD2.mat
!rm ../../inputs/WGTSMACE/WGTS_CDProfileIntegral/LambdaIntegral_Table-WGTSProfile_100-WGTSzCells_1.1e+17MolPerCm2.mat
[SamakTe SamakRFfpd] = BuildfpdSamakRF();
end

%% Definition
E0    = 18575;
qu    = [18373  18483  18488  18493  18498  18503  18508  18513  18518  18523  18528  18533  18535  18537  18539  18541  18543  18545  18547  18549  18551  18553  18555  18557  18559  18561  18562  18563  18564  18565  18566  18567  18569  18571  18573  18578  18583  18593  18603  18623];
qUbin = 10;

close all

%% SuperImpose Master / Samak
fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
r=plot(MasterTe,MasterRFfpd(:,qUbin),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
hold on
s=plot(SamakTe,SamakRFfpd(:,qUbin),'-','Color',rgb('IndianRed'),'LineWidth',2);
hold off
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e(eV)'),'FontSize',16);
xtickangle(45);xtickformat('%0.4f');
ylabel('transmission','FontSize',16);
legend([r s],'Master RF','Samak RF','Location','NorthWest');
xlim([qu(qUbin)-5 SamakTe(end)]);
grid on
PrettyFigureFormat
set(gca,'FontSize',18);
export_fig(gcf,'KNM1_Samak_Master_RFcomparison_GoldenPixelList_1.png','-m3');

% SuperImpose Master / Samak - Interpolated
fMasterRFfpd = @(te) interp1(MasterTe,MasterRFfpd(:,qUbin),te);
fSamakRFfpd = @(te) interp1(SamakTe,SamakRFfpd(:,qUbin),te);
te=18320:0.1:18600;

fign = figure('Name','Samak','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
maintitle=sprintf('Master Response Function - \\rho.d = %.4g mol/cm^2 - Bs = %.2f T - Ba = %.3g G - qU=-%.2f V',...
    1.1e17,2.52,6.3,qu(qUbin));
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';

s1=subplot(2,1,1)

r=stairs(te,fMasterRFfpd(te),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
hold on
s=stairs(te,fSamakRFfpd(te),'-','Color',rgb('IndianRed'),'LineWidth',2);
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e (eV)'),'FontSize',18);
ylabel('transmission','FontSize',18);
legend([r s],'Master RF - Interpolated','Samak RF - Interpolated','Location','NorthWest');
grid on
PrettyFigureFormat
xlim([qu(qUbin) qu(qUbin)+3]);
set(gca,'FontSize',18);

s2=subplot(2,1,2)
r=stairs(te,(fSamakRFfpd(te)-fMasterRFfpd(te)),'-','Color',rgb('DarkSlateGray'),'LineWidth',2);
set(gca, 'YScale', 'lin');
xlabel(sprintf('T_e  (eV)'),'FontSize',18);
xlim([qu(qUbin) qu(qUbin)+3]);
ylabel('transmission','FontSize',18);
legend([r],'Samak - Master RF','Location','NorthWest');
grid on
PrettyFigureFormat
set(gca,'FontSize',18);
linkaxes([s1,s2],'x');
export_fig(gcf,'KNM1_Samak_Master_RFcomparison_GoldenPixelList_2.png','-m3');
% Plot Sensitivity for 3 ranges
% as a function of Ba (and Background)
% loop over systematics
clear;
MACE_Ba_T_all = (3:12).*1e-04;         % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
range_all = [30,45,60];    % MTD energy range: 30,45,60eV below E0 (18575eV)
WGTS_B_T = 0.7*3.6;
AnchorBkg6G = 0.238;
SysBudget = '03';  % Systematics Input Information in Result folder
FPDeff = 0.9;
BKG = 0.210;
PlotFontSize = 18;
nSys = 5;
TimeSec = 3*365*24*60*60;
%Init: gather results
mNu90     = zeros(numel(range_all),numel(MACE_Ba_T_all),nSys);                  % sensitivity on mnu^2 (90% C.L.)
mNumin    = zeros(numel(range_all),numel(MACE_Ba_T_all),nSys);                  % mass with minimal chi2

InfoBox = 'OFF'; % title, annotations

for i=1:numel(range_all)
    range = range_all(i);
    for b=1:numel(MACE_Ba_T_all)
        MACE_Ba_T = MACE_Ba_T_all(b);
        TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
        save_name = sprintf('./results/%s_%0.0fd_SensitivityNominal_ResultsNuMassScan_MTD-%s_Bs%.2fT_FPDeff%.2f_BKG-%.0fmcps_Systematics.mat',SysBudget,TimeSec/86400,strrep(TD,'Sensitivity_',''),WGTS_B_T,FPDeff,BKG*1e3);
        if exist(save_name,'file')==2
            d = importdata(save_name);
        else
            fprintf('File doesnt exist! \n');
            return
        end
        
        mNu90(i,b,:)  = d.mNu90(1:nSys); %last one is all systematics combined
        mNumin(i,b,:) = d.mNumin(1:nSys);
    end
end
SysEffects = d.mySysEffects;
%% Plot overview
mNu90max = 1e3*max(max(max(sqrt(mNu90))));
mNu90min = 1e3*min(min(min(sqrt(mNu90)))); %200;

BKG = GetBackground('MACE_Ba_T',MACE_Ba_T_all,'WGTS_B_T',WGTS_B_T,'AnchorBkg6G',AnchorBkg6G);
f80 = figure('Name','Overview','Renderer','opengl','Visible','off');
maintitle = ['m_{\nu} ',sprintf('Sensitivity (Scan Method) for KATRIN nominal settings (%.0f %% B_s and B_{max})',WGTS_B_T/3.6*100)];
set(f80, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
if strcmp(InfoBox,'ON')
a=annotation('textbox', [0 0.9 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=12;a.FontWeight='bold';
end
s1 = subplot(3,2,[1.48,1.57]);
SysIndex = 1;
pnone = plot(MACE_Ba_T_all.*1e4,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
hold on;
p30 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(1,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('CornFlowerBlue'));
p45 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(2,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('ForestGreen'));
p60 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(3,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('IndianRed'));
legend([pnone,p30,p45,p60],'stat','E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
legend boxoff
%xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
xlim([min(MACE_Ba_T_all.*1e4) max(MACE_Ba_T_all.*1e4)])
PrettyFigureFormat;
set(gca,'FontSize',14);
ylabel('L({m_{\nu}}) 90% C.L. (meV)','FontSize',13);
grid on;
ylim([mNu90min mNu90max]);
%Get Background axes
axBKG = gca;
ax1_pos = axBKG.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlim([min(MACE_Ba_T_all.*1e4) max(MACE_Ba_T_all.*1e4)]);
xlabel('Background (mcps)');
xticks(MACE_Ba_T_all.*1e4);
yticks([]);
xticklabels(string(round(BKG*1e3)));
PrettyFigureFormat;
set(gca,'FontSize',14);

%%
close all
f80 = figure('Name','Overview','Renderer','opengl');
set(f80, 'Units', 'normalized', 'Position', [0, 0, 1.5, 1.5]);
for i=1:nSys-1 %with systematics
    switch SysBudget
        case {'01','03'}
            s2 = subplot_er(2,2,i);%+2
      
        case '00'
            s2 = subplot_er(2,2,i+1);
    end
    SysIndex = i+1;
    pnone = plot(MACE_Ba_T_all.*1e4,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
    none_leg = sprintf('stat + sys (%s)',SysEffects{SysIndex-1});
    hold on;
    p30 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(1,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('CornFlowerBlue'));
    p45 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(2,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('ForestGreen'));
    p60 = plot(MACE_Ba_T_all.*1e4,1e3.*sqrt(mNu90(3,:,SysIndex)),'o--','LineWidth',3,'Color',rgb('IndianRed'));
    %legend([pnone,p30,p45,p60],none_leg,'E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
    leg = legend(pnone,none_leg,'Location','best');
    legend boxoff
    xlim([min(MACE_Ba_T_all.*1e4) max(MACE_Ba_T_all.*1e4)])
    PrettyFigureFormat; grid on;
    set(gca,'FontSize',PlotFontSize);
        if SysIndex==2
        ylabel('L(m_{\nu}) 90% C.L. (meV)','FontSize',PlotFontSize);
    elseif SysIndex==3
    elseif SysIndex==4
        ylabel('L(m_{\nu}) 90% C.L. (meV)','FontSize',PlotFontSize);
        xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
    elseif SysIndex==5
        xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
        end
    leg.FontSize = PlotFontSize;
    ylim([mNu90min mNu90max]);
  
end
linkaxes([s1,s2],'xy');
set(gcf,'Visible','on');
if ~exist('../sensitivity_nominalKATRIN/plots/png/BaDependency','dir')
    mkdir ../sensitivity_nominalKATRIN/plots/png/BaDependency/
    mkdir ../sensitivity_nominalKATRIN/plots/pdf/BaDependency/
    mkdir ../sensitivity_nominalKATRIN/plots/fig/BaDependency/
end
save_name = sprintf('%s_SensitivityNominal_NuMassScan_BKG_Ba_%.0fBFields_Overview',SysBudget,WGTS_B_T/3.6*100);
print(f80,['../sensitivity_nominalKATRIN/plots/png/BaDependency/',save_name,'.png'],'-dpng');
savefig(f80,['../sensitivity_nominalKATRIN/plots/fig/BaDependency/',save_name,'.fig'],'compact');
publish_figurePDF(f80,['../sensitivity_nominalKATRIN/plots/pdf/BaDependency/',save_name,'.pdf']);


%% stat only + sensitivity of TDR
close
f81 = figure('Name','TDR','Renderer','opengl');
set(f81, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
x = MACE_Ba_T_all*1e4;
pnone = plot(x,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
hold on;
p30 = plot(x,1e3.*sqrt(sqrt(mNu90(1,:,1).^2+(1.64*0.017).^2)),'o--','LineWidth',3,'Color',rgb('CornFlowerBlue'));
p45 = plot(x,1e3.*sqrt(sqrt(mNu90(2,:,1).^2+(1.64*0.017).^2)),'o--','LineWidth',3,'Color',rgb('ForestGreen'));
p60 = plot(x,1e3.*sqrt(sqrt(mNu90(3,:,1).^2+(1.64*0.017).^2)),'o--','LineWidth',3,'Color',rgb('IndianRed'));
xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
%xlabel('Background (mcps)')
ylabel('\sigma_{m_{\nu}} 90% C.L. (meV)');
xlim([min(x) max(x)]);
ylim([150 375]);
yticks([150:25:375]);
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);
grid on;
legend([pnone,p30,p45,p60],'stat + sys (\sigma=0.017eV^2)','E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
legend boxoff
ylim([mNu90min mNu90max]);

%Get Background axes
axBKG = gca;
ax1_pos = axBKG.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlabel('Background (mcps)');
xlim([min(x) max(x)]);
xticks(x);
yticks([]);
xticklabels(string(round(BKG*1e3)))
PrettyFigureFormat
set(gca,'FontSize',PlotFontSize);
% 
% save_name = sprintf('%s_SensitivityNominal_NuMassScan_BKG_Ba_%.0fBFields_TDRsystematics',SysBudget,WGTS_B_T/3.6*100);
% export_fig(f81,['../sensitivity_nominalKATRIN/plots/png/BaDependency/',save_name,'.png']);
% savefig(f81,['../sensitivity_nominalKATRIN/plots/fig/BaDependency/',save_name,'.fig'],'compact');
% publish_figurePDF(f81,['../sensitivity_nominalKATRIN/plots/pdf/BaDependency/',save_name,'.pdf']);
% 

%% plot 1 effect
close
InfoBox = 'OFF';
ThisFontSize = 26;
SysIndex = 1; %1=stat, 2=TC, 3=FSD, 4=RF,5=all
f82 = figure('Name','SysEffects','Renderer','opengl');
%set(f82, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.8]);
set(f82, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.8]);
x = MACE_Ba_T_all*1e4;
pnone = plot(x,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
hold on;
p30 = plot(x,1e3.*sqrt(mNu90(1,:,SysIndex)),'o--','LineWidth',7,'Color',rgb('CornFlowerBlue'));
p45 = plot(x,1e3.*sqrt(mNu90(2,:,SysIndex)),'o--','LineWidth',7,'Color',rgb('ForestGreen'));
p60 = plot(x,1e3.*sqrt(mNu90(3,:,SysIndex)),'o--','LineWidth',7,'Color',rgb('IndianRed'));
xlabel('B-field in analyzing plane B_a (10^{-4} T)')
ylabel('L(m_{\nu}) 90% C.L. (meV)');
xlim([min(MACE_Ba_T_all.*1e4) max(MACE_Ba_T_all.*1e4)])
ylim([mNu90min mNu90max]);
yticks([250:50:400]);
PrettyFigureFormat;
set(gca,'FontSize',ThisFontSize);
grid on;
ylim([mNu90min mNu90max]);
if SysIndex==1
none_leg = 'stat';
SysName = 'stat';
else
none_leg = sprintf('stat + sys (%s)',SysEffects{SysIndex-1});
SysName = SysEffects{SysIndex-1};
end
legend([pnone,p30,p45,p60],none_leg,'E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
legend boxoff
leg.FontSize = ThisFontSize;
axBKG = gca;
ax1_pos = axBKG.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlabel(' ');
xlim([min(x) max(x)]);
xticks(x);
yticks([]);
xticklabels(string(round(BKG*1e3)))

FontName = 'Helvetica';
FontSize = ThisFontSize;
FontWeight = 'bold';
set(gca,'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight)
set(get(gca,'XLabel'),'FontName',FontName,'FontSize',FontSize,'FontWeight',FontWeight);
set(gca,'FontSize',ThisFontSize);
set(gca,'TickLength',[.01 .01]);
set(gca,'XMinorTick','on');

if strcmp(InfoBox,'ON')
if strcmp(SysBudget,'03')
   MACE_Ba_Err=sprintf('(B_a \\cdot 3.9\\cdot10^{-3}+1.07\\cdot 10^{-6})'); %absolute
   Ba_unit=sprintf('');
else
    MACE_Ba_Err = 1e6*s.MACE_Ba_T_Err;
    Ba_unit = sprintf('\\muT');
end

if SysIndex==3 %FSD
    SysUncertainties =  [sprintf('\nFinal State Distribution:'),...
        sprintf(' Normalization %.0f %%, ',s.FSDNorm_RelErr*100),...
        sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*s.FSDShapeGS_RelErr,100*s.FSDShapeES_RelErr)];
elseif SysIndex==4 %RF
    SysUncertainties = [sprintf('Response Function: '),...
        sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
        sprintf('\\DeltaB_a = %s , \\DeltaB_{T/max} = %.1f%%, ',s.MACE_Ba_T_Err_str,100*s.MACE_Bmax_T_RelErr),...
        sprintf('\nEnergy Loss Uncertainties from Eur. Phys. J. D 10, 39–52')];
elseif SysIndex==5 %all
    SysUncertainties = [sprintf('Response Function: '),...
        sprintf('\\Delta \\rhod\\sigma = %.1f%%, ',100*sqrt(s.WGTS_CD_MolPerCm2_RelErr^2+s.ISXsection_RelErr^2)),...
        sprintf('\\DeltaB_a = %s %s, \\DeltaB_{T/max} = %.1f%%, ',MACE_Ba_Err,Ba_unit,100*s.MACE_Bmax_T_RelErr),...
        sprintf('\nEnergy Loss Uncertainties from Eur. Phys. J. D 10, 39–52'),...
        sprintf('\nFinal State Distribution:'),...
        sprintf(' Normalization %.0f %%, ',s.FSDNorm_RelErr*100),...
        sprintf('Bin-to-Bin uncorrelated %.0f %% (GS), %.0f %% (ES)',100*s.FSDShapeGS_RelErr,100*s.FSDShapeES_RelErr)];
else
    SysUncertainties = '';
end

a=annotation('textbox', [0.14 0.12 1 0.1], ...
    'String', SysUncertainties, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left');
a.FontSize=11;a.FontWeight='bold';
end
maintitle = 'background (mcps)';
a=annotation('textbox', [0 0.9 1 0.02], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=ThisFontSize;a.FontWeight='bold';


save_name = sprintf('%s_SensitivityNominal_NuMassScan_BKG_Ba_%.0fBFields_%s',SysBudget,WGTS_B_T/3.6*100,SysName);
export_fig(f82,['../sensitivity_nominalKATRIN/plots/png/BaDependency/',save_name,'.png']);
savefig(f82,['../sensitivity_nominalKATRIN/plots/fig/BaDependency/',save_name,'.fig'],'compact');
publish_figurePDF(f82,['../sensitivity_nominalKATRIN/plots/pdf/BaDependency/',save_name,'.pdf']);

%% plot sys RF
close
SysIndex = 4;%4=RF
f83 = figure('Name','RFSys','Renderer','opengl');
set(f83, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
x = MACE_Ba_T_all*1e4;
SysErr30 = sqrt(sqrt(mNu90(1,:,SysIndex).^2-mNu90(1,:,1).^2))*1e3;
SysErr45 = sqrt(sqrt(mNu90(2,:,SysIndex).^2-mNu90(1,:,1).^2))*1e3;
SysErr60 = sqrt(sqrt(mNu90(3,:,SysIndex).^2-mNu90(1,:,1).^2))*1e3;

pnone = plot(x,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
hold on;
p30 = plot(x,SysErr30,'o--','LineWidth',3,'Color',rgb('CornFlowerBlue'));
p45 = plot(x,SysErr45,'o--','LineWidth',3,'Color',rgb('ForestGreen'));
p60 = plot(x,SysErr60,'o--','LineWidth',3,'Color',rgb('IndianRed'));
xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
ylabel('\sigma_{sys_{m_{\nu}}} 90% C.L. (meV)');
xlim([min(MACE_Ba_T_all.*1e4) max(MACE_Ba_T_all.*1e4)])
%ylim([150 500]);
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);
grid on;
none_leg = sprintf('\\sigma_{sys(%s)}^2=\\sigma_{stat+sys}^2-\\sigma_{stat}^2 ',SysEffects{SysIndex-1});
leg = legend([p30,p45,p60],'E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
leg.Title.String = none_leg;
legend boxoff
axBKG = gca;
ax1_pos = axBKG.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlabel('Background (mcps)');
xlim([min(x) max(x)]);
xticks(x);
yticks([]);
xticklabels(string(round(BKG*1e3)))
PrettyFigureFormat
set(gca,'FontSize',PlotFontSize);

save_name = sprintf('%s_SensitivityNominal_NuMassScan_BKG_Ba_%.0fBFields_%s_RFsystematic',SysBudget,WGTS_B_T/3.6*100,SysEffects{SysIndex-1});
export_fig(f83,['../sensitivity_nominalKATRIN/plots/png/BaDependency/',save_name,'.png']);
savefig(f83,['../sensitivity_nominalKATRIN/plots/fig/BaDependency/',save_name,'.fig'],'compact');
publish_figurePDF(f83,['../sensitivity_nominalKATRIN/plots/pdf/BaDependency/',save_name,'.pdf']);

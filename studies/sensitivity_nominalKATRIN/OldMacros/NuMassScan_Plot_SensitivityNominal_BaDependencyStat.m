% Plot Sensitivity for 3 ranges
% as a function of Ba (and Background)
% loop over systematics
clear;
MACE_Ba_T_all = (3:12).*1e-04;    % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
WGTS_B_T = 3.6*0.7;
range_all = [30,45,60];    % MTD energy range: 30,45,60eV below E0 (18575eV)
%Init: gather results
mNu90     = zeros(numel(range_all),numel(MACE_Ba_T_all));                  % sensitivity on mnu^2 (90% C.L.)
mNumin    = zeros(numel(range_all),numel(MACE_Ba_T_all));                  % mass with minimal chi2

for i=1:numel(range_all)
    range = range_all(i);
    for b=1:numel(MACE_Ba_T_all)
        MACE_Ba_T = MACE_Ba_T_all(b);
        TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
        save_name = sprintf('./results/SensitivityNominal_ResultsNuMassScan_MTD%s_Bs%.2fT_Stat.mat',strrep(TD,'Sensitivity_',''),WGTS_B_T);
        if exist(save_name,'file')==2
            d = importdata(save_name);
        else
            fprintf('File doesnt exist! \n');
            return
        end
        
        mNu90(i,b)  = d.mNu90; %last one is all systematics combined
        mNumin(i,b) = d.mNumin;
    end
end

%% Plot
close
x = BKG*1e3;
%x=MACE_Ba_T_all.*1e4;
BKG = GetBackground('MACE_Ba_T',MACE_Ba_T_all,'WGTS_B_T',WGTS_B_T);
f80 = figure(80);
maintitle = ['m_{\nu} ',sprintf('Sensitivity (Scan Method) for KATRIN nominal settings')];
set(f80, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
  a=annotation('textbox', [0 0.81 0.73 0.1], ...
            'String', maintitle, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        a.FontSize=14;a.FontWeight='bold';
pnone = plot(x,NaN.*zeros(numel(MACE_Ba_T_all.*1e4),1),'Color',rgb('White'));
hold on;
p30 = plot(x,sqrt(sqrt(mNu90(1,:).^2+0.012^2)),'o--','LineWidth',3,'Color',rgb('CadetBlue'));
p45 = plot(x,sqrt(sqrt(mNu90(2,:).^2+0.012^2)),'o--','LineWidth',3,'Color',rgb('IndianRed'));
p60 = plot(x,sqrt(sqrt(mNu90(3,:).^2+0.012^2)),'o--','LineWidth',3,'Color',rgb('DarkGoldenrod'));
legend([pnone,p30,p45,p60],'stat','E_{0}-30 eV','E_{0}-45 eV','E_{0}-60 eV','Location','best');
legend boxoff
ylabel('\sigma_{m_{\nu}} 90% C.L. (eV)');
if all(x==MACE_Ba_T_all*1e4)
    xlabel('B-Field in Analyzing Plane B_a (10^{-4} T)')
elseif all(x==BKG*1e3)
    xlabel('Background (mcps)');
end

xlim([min(x) max(x)])
PrettyFigureFormat;
set(gca,'FontSize',16);
grid on;
%Background Axes
axBKG = gca;
%ax1_pos = axBKG.Position; % position of first axes
%ax2 = axes('Position',ax1_pos,...
%    'XAxisLocation','top',...
 %   'YAxisLocation','right',...
%    'Color','none');
%xlabel('Background (mcps)');
%hold on;
%plot(BKG*1e3,NaN*zeros(numel(BKG),1),'parent',ax2);
%xlim([min(BKG*1e3) max(BKG*1e3)]);
%yticklabels = '';
%ax2.YTick = [];
%set(ax2,'Xdir','reverse');
%PrettyFigureFormat; 
%set(gca,'FontSize',16);



save_name = sprintf('SensitivityNominal_NuMassScan_BKG_Ba');
export_fig(f80,['./plots/png/',save_name,'.png']);
savefig(f80,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f80,['./plots/pdf/',save_name,'.pdf']);
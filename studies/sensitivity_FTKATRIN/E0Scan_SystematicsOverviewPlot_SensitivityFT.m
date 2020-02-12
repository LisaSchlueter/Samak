% Overview plot for 1setting - all covariance matrices
% Settings
mySysEffects  = {'TC','TASR','FSD','RF','all'};
Q_i_Fit       = -0.75:0.01:0.75; % good for  400eV
%Q_i_Fit       = -0.95:0.01:0.95; % good for  200eV
%Q_i_Fit        = -0.4:0.01:0.4; % good for 1600eV
nQi           = numel(Q_i_Fit);
exclDataStart = 7;
RunList       = 'StackCD100all';
Time          = sprintf('5.2 days'); % for StackCD100all
RecomputeFlag = 'ON';                % Compute everything again - do not load form file

switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7
        belowE0 = 400;
    case 9
        belowE0 = 200;
end
save_name = sprintf('./results/SensitivityFT_ResultsE0Scan_%s_%.0feVrange.mat',RunList,belowE0);
if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
    load(save_name);
else
%Init: gather results
chi2min   = zeros(nQi,numel(mySysEffects)+1);
E0min     = zeros(numel(mySysEffects)+1,1);
E090low   = zeros(numel(mySysEffects)+1,1);
E090up    = zeros(numel(mySysEffects)+1,1);

[~, chi2min(:,1), dof, E0min(1), E090low(1), E090up(1)] = E0Scan_SensitivityStudy_FTKATRIN(...
    'chi2','chi2Stat','Q_i_Fit',Q_i_Fit,'exclDataStart',exclDataStart,'RunList',RunList);
parfor i=1:numel(mySysEffects)
 [~, chi2min(:,i+1),~, E0min(i+1), E090low(i+1), E090up(i+1)] = E0Scan_SensitivityStudy_FTKATRIN(...
    'chi2','chi2CM','SysEffect',mySysEffects{i},'Q_i_Fit',Q_i_Fit,'exclDataStart',exclDataStart,'RunList',RunList);   
end

if exist('./results/','dir')~=7 %in case results folder doesnt exist
    mkdir ./results/
end
save(save_name,'chi2min','dof','E090low','E090up','E0min','Q_i_Fit');
end
%% Plot
close 
% extract systematic uncertainty
SysErrlow = sqrt(E090low.^2-E090low(1)^2)*1e3;
SysErrup  = sqrt(E090up.^2-E090up(1)^2)*1e3;
SysErrTC_str = ['\sigma_{sysTC} = ',sprintf('- %.0f + %.0f meV',SysErrlow(2),SysErrup(2))];
SysErrTASR_str = ['\sigma_{sysTASR} = ',sprintf('- %.0f + %.0f meV',SysErrlow(3),SysErrup(3))];
SysErrFSD_str = ['\sigma_{sysFSD} = ',sprintf('- %.0f + %.0f meV',SysErrlow(4),SysErrup(4))];
SysErrRF_str = ['\sigma_{sysRF} = ',sprintf('- %.0f + %.0f meV',SysErrlow(5),SysErrup(5))];
SysErrall_str = ['\sigma_{sysAll} = ',sprintf('- %.0f + %.0f meV',SysErrlow(6),SysErrup(6))];
SysErr_str = {SysErrTC_str,SysErrTASR_str,SysErrFSD_str,SysErrRF_str,SysErrall_str};

%legend entries
stat_leg = ['Stat:                                                         \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(1))*1e3, E090up(1)*1e3)];
TC_leg = ['Stat + Theoretical Corrections CM:        \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(2))*1e3, E090up(2)*1e3)];
TASR_leg = ['Stat + Tritium Activity Fluctuatoin CM:  \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(3))*1e3, E090up(3)*1e3)];
FSD_leg = ['Stat + Final State Distribution CM:        \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(4))*1e3, E090up(4)*1e3)];
RF_leg = ['Stat + Response Function CM:                \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(5))*1e3, E090up(5)*1e3)];
All_leg = ['Stat + Combined CM:                              \sigma_{E0eff} = ',sprintf('- %.0f + %.0f meV',abs(E090low(6))*1e3, E090up(6)*1e3)];

f33 = figure(3);
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
pScan = plot(Q_i_Fit,chi2min,'LineWidth',3);
hold on;
ptitle = plot(Q_i_Fit,NaN*zeros(nQi,1),'Color',rgb('White'));
leg = legend([ptitle; pScan],'Sensitivity on effective Endpoint (90% C.L.)',stat_leg,TC_leg,TASR_leg,FSD_leg,RF_leg,All_leg);
leg.NumColumns = 1;
leg.Location = 'north';
legend boxoff
xlabel('18573.7eV - E_{0eff} (eV)')
ylabel(['\chi^2 (',num2str(dof),' dof)']);
xlim([min(Q_i_Fit) max(Q_i_Fit)]);
ylim([0 max(max(chi2min))]);
PrettyFigureFormat;
set(gca,'FontSize',18);
leg.FontSize = 14;
title(sprintf('Samak Fit to Simulation (asimov) \n - KATRIN First Tritium %s - MTD: %s  - %.0feV range',...
   Time,RunList,belowE0));
a = annotation('textbox',[0.4 0.25 0.2 0.4],'String',SysErr_str,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=14; a.FontWeight='bold';
%% save plot
if ~exist('./plots/png','dir')
    mkdir ./plots/png/
    mkdir ./plots/pdf/
    mkdir ./plots/fig/
end
savePlot_name = sprintf('SensitivityFT_E0ScanOverview_%s_%.0feV',RunList,belowE0);
print(f33,['./plots/png/',savePlot_name,'.png'],'-dpng');
savefig(f33,['./plots/fig/',savePlot_name,'.fig'],'compact');
publish_figurePDF(f33,['./plots/pdf/',savePlot_name,'.pdf']);
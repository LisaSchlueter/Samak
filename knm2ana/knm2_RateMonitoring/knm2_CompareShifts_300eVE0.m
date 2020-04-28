% comparison plot for "shifts"
% 300eV rate analysis versus large range E0 fit

Baseline = ' ' ; %'RW2PSR2';

%% load 300 eV analysis
file300eV = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/knm2_RManalysis_Ring.mat'];

if exist(file300eV,'file')
    d = importdata(file300eV);
    Shift300eV    = d.Shift;
    Shift300eVErr = d.ShiftErr;
else
    ROIstr = 'Default';
    [~, ~, Shift300eV, Shift300eVErr] = knm2_RingwiseFit('ROI',ROIstr);
end
%% load large range E0 fit
[Diffs_E0, E0, E0Err,FitResultsLargeRange,~,~,~] = knm2_LargeRangeFit_Ringwise;

if strcmp(Baseline,'Mean')
    Diffs_E0 = Diffs_E0-mean(mean(Diffs_E0));
    Shift300eV = Shift300eV - mean(mean(Shift300eV));
else
    Baseline = 'RW2PR1';
end
% Shift300eV = Shift300eV-mean(mean(Shift300eV));
% Diffs_E0 = Diffs_E0-mean(mean(Diffs_E0));
%Shift300eVErr = Shift300eVErr.*sqrt([171;97;93]);%sqrt(361);
%% plot result
fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[0.1, 0.1,0.5,0.8]);
PlotStyle = { 'o','MarkerSize',6,'LineWidth',2,'CapSize',0};

ymin1 = min(min(-Diffs_E0-E0Err)).*1e3;
ymax1 = max(max(-Diffs_E0+E0Err)).*1e3;
ymin2 = min(min(Shift300eV-Shift300eVErr)).*1e3;
ymax2 = max(max(Shift300eV+Shift300eVErr)).*1e3;
ymin = 1.05*((ymin1<ymin2).*ymin1+(ymin1>ymin2).*ymin2);
ymax = 1.8*((ymax1>ymax2).*ymax1+(ymax1>ymax2).*ymax2);

for i = 1:3
    subplot(3,1,i);
    e1 = errorbar(-Diffs_E0(i,:).*1e3,E0Err(i,:).*1e3,PlotStyle{:},...
        'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
    hold on;
    e2 = errorbar(Shift300eV(i,:).*1e3,Shift300eVErr(i,:).*1e3,PlotStyle{:},...
        'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
    xlabel(sprintf('Pseudo-ring'));
    ylabel(sprintf('U - \\langleU\\rangle (mV)'));
    leg = legend(sprintf('{\\itE}_0 fit [-90 to -45, +135] eV'),sprintf('300 eV rate analysis'),...
        'Location','northwest','EdgeColor',rgb('Silver'));
    leg.Title.String = sprintf('RW period %.0f',i);
    leg.Title.FontWeight = 'normal';
    xlim([0.5 4.5]);
    
   % ylim([ymin, ymax]);
    
    PrettyFigureFormat('FontSize',20);
    xticks([1:4]);
    set(gca,'XMinorTick','off');
end

plotdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
plotname = sprintf('%sknm2_CompareShifts_300E0_Basline%s.png',plotdir,Baseline);
if contains(plotname,'.pdf')
export_fig(fig88,plotname);
elseif contains(plotname,'.png')
    print(fig88,plotname,'-dpng','-r500');
end
fprintf('save plot to %s \n',plotname)


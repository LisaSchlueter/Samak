% Loop over Sys Effects
% Comarison of 2 methods + Fits to Data
RunList = 'StackCD100all';
nSamples = 5000;
fixPar = '1 5 6';
Q_i    = 18573.7;
exclDataStart = 7;
mySysEffects  = {'TC','TASR','FSD','RF','all'};
switch exclDataStart
    case 7
        belowE0 = 402;
    case 9
        belowE0 = 202;
end

%% Load MC sensitivities
E090low_MC = zeros(numel(mySysEffects)+1,1);
E090up_MC = zeros(numel(mySysEffects)+1,1);
%stat only
file_nameMC = sprintf('./results/SensitivityStudy_FTKATRIN_%s_StatFluctON_SysFluctOFF_chi2Stat_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
    RunList,strrep(fixPar,' ',''),belowE0,Q_i*10,nSamples);
%result_MC = importdata(file_nameMC);
% E0      = sort(cellfun(@(x) x.par(2),result_MC.FitResult)+Q_i);
% E090low_MC(1) = E0(nSamples*0.05)-mean(E0); % Sensitivity 90% C.L.
% E090up_MC(1)  = E0(nSamples*0.95)-mean(E0);


% %stat + sys
% for i=1:numel(mySysEffects)
%     SysEffect = mySysEffects{i};
%     file_nameMC = sprintf('./results/MC_SensitivityStudy_FTKATRIN_%s_StatFluctON_SysFluctON_chi2CM%s_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
%         RunList,SysEffect,strrep(fixPar,' ',''),belowE0,Q_i*10,nSamples);
%     result_MC = importdata(file_nameMC);
%     E0      = sort(cellfun(@(x) x.par(2),result_MC.FitResult)+Q_i);
%     E090low_MC(i+1) = E0(nSamples*0.05)-mean(E0); % Sensitivity 90% C.L.
%     E090up_MC(i+1)  = E0(nSamples*0.95)-mean(E0);
% end

%% Load Scan sensitivies
file_nameScan = sprintf('./results/SensitivityFT_ResultsE0Scan_%s_%.0feVrange.mat',RunList,belowE0-2);
result_Scan = importdata(file_nameScan);
E090low_Scan = result_Scan.E090low;
E090up_Scan = result_Scan.E090up;

%% Load Uncertainty from Real Data: attention: here E0err is 1 sigma C.L.
file_nameFit = sprintf('./../ft_Stacking/results/ResultE0_SystematicBreakdown_%s_%.0feV.mat',RunList,belowE0);
result_Fit = importdata(file_nameFit);
E090_Fit = result_Fit.E0Err.*1.64;
%% Plot
close;
f19 = figure('Renderer','opengl');
set(f19, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
x=(1.0:1:6.0);
%plot Scan
pScan = plot(x,abs(E090up_Scan),'o--','LineWidth',3,'MarkerSize',10,...
    'Color',rgb('CadetBlue'),'MarkerFaceColor',rgb('CadetBlue'));
hold on;
plot(x,E090low_Scan,'o--','LineWidth',3,'MarkerSize',10,...
    'Color',rgb('CadetBlue'),'MarkerFaceColor',rgb('CadetBlue'));
plot(x,zeros(numel(x),1),'k-','LineWidth',2);
%plot MC
%pMC = plot(x,E090low_MC,'o--','LineWidth',2,'MarkerSize',10,...
 %   'Color',rgb('Goldenrod'),'MarkerFaceColor',rgb('GoldenRod'));
%plot(x,E090up_MC,'o--','LineWidth',2,'MarkerSize',10,...
%    'Color',rgb('Goldenrod'),'MarkerFaceColor',rgb('GoldenRod'));
%plot Fit
pFit = plot(x,E090_Fit,'o--','LineWidth',3,'MarkerSize',10,...
    'Color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'));
xticklabels({'Stat.','Theoretical Corr.','Tritium Activity Fluct.','FSD','Response Function','Combined'})
ylabel(sprintf('sensitivity E_{0eff} 90%% C.L. (eV)'))
grid on;
PrettyFigureFormat
ax = gca;
ax.XTick = x; %pMC   % sprintf('Monte Carlo (%.0f Samples)',nSamples)
legend([pScan,pFit],'Asimov Scan Method (\chi^2+2.7)',...
    'Uniform Fit to Data','Location','northwest');
legend boxoff;
set(gca,'FontSize',18);
ylim([0 0.45])
title(sprintf('Sensitivity on effective Endpoint - KATRIN First Tritium \n MTD: %s (5.2 days) - %.0feV scan range - \\nu-mass fixed',RunList, belowE0-2));

save_name = sprintf('E0SensitivityOverview_%s_%0.feVrange',RunList,belowE0-2);
export_fig(f19,['./plots/png/',save_name,'.png']);
savefig(f19,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f19,['./plots/pdf/',save_name,'.pdf']);

%% BarPlot
f333 = figure('Renderer','opengl');
set(f333, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.7]);
pScan = bar(x-0.15,abs(E090up_Scan),'FaceColor',rgb('CadetBlue'),'BarWidth',0.3);
pScan.LineStyle = 'none';
hold on;
pFit = bar(x+0.15,E090_Fit,'FaceColor',rgb('FireBrick'),'BarWidth',0.3);
pFit.LineStyle = 'none';
grid on;
PrettyFigureFormat
set(gca,'FontSize',14);
ax = gca;
ax.XTick = x;
leg = legend([pScan,pFit],'Simulation (\chi^2+2.7)',...
    'Uniform Fit to FT Data','Location','northwest');
legend boxoff;
leg.FontSize = 18;
xticklabels({'Stat.','Theoretical Corr.','Tritium Activity Fluct.','FSD','Response Function','Combined'})
ylabel(sprintf('sensitivity E_{0eff} 90%% C.L. (eV)'),'FontSize',18)
ylim([0 0.45])
xlim([0.5 6.5]);

title(sprintf('Sensitivity on effective Endpoint - KATRIN First Tritium \n MTD: %s (5.2 days) - %.0feV scan range - \\nu-mass fixed',RunList, belowE0-2));

save_name = sprintf('E0SensitivityOverview_Bar_%s_%0.feVrange',RunList,belowE0-2);
export_fig(f333,['./plots/png/',save_name,'.png']);
savefig(f333,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f333,['./plots/pdf/',save_name,'.pdf']);

% Plot default bin sizes of FSDs
Arg = {'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'DataType','Real',...
    'RunList','KNM2_Prompt'};

% get run-wise endpoint
D = MultiRunAnalysis(Arg{:});
exE     = D.ModelObj.TTexE;
binSize = diff(exE(2:end-1));

%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
b1 = bar(exE(2:end-2),binSize,'BarWidth',7,'FaceColor',rgb('DodgerBlue'));
xlabel('Excitation energy (eV)');
ylabel('Bin size (eV)');
PrettyFigureFormat('FontSize',24);
xlim([0,40]);
ylim([0 3.1])
leg = legend(b1,sprintf('T_2 FSD: %s',strrep(D.FSDFlag,'K',' K')));
legend boxoff;
leg.FontSize = get(gca,'FontSize');
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/plots/'];
savename = [savedir,sprintf('FSDbinning_%s.pdf',D.FSDFlag)];
export_fig(f1,savename);



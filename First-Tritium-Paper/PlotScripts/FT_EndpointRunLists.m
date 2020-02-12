% endpoint with different run lists
%% set up models
Mode = 'RD';
exclDataStart = 13;
chi2 = 'chi2Stat';
RecomputeFlag = 'OFF';
switch Mode
    case 'Scan'
        RunList1 = 'StackCD100down';
        RunList2 = 'StackCD100up';
        RunList3 = 'StackCD100random';
        RunList4 = 'FTpaper';
        myXlabel = 'Scan direction';
        myXticks = {'Down','Up','Random','All'};    
    case 'RD'
        RunList1 = 'FT_RD24';
        RunList2 = 'FT_RD48';
        RunList3 = 'FT_RD72';
        RunList4 = 'FTpaper';
        myXlabel = 'Column density (%)';
     case 'Scan2'
        RunList1 = 'StackAnnaDown';
        RunList2 = 'StackAnnaUp';
        RunList3 = 'StackAnnaRandom';
        RunList4 = 'FTpaper';
        myXlabel = 'Scan direction';
        myXticks = {'Down','Up','Random','All'};           
end

savedir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/'];
savename = [savedir,sprintf('results/FT_EndpointRunLists_%.0f_%s_%s.mat',exclDataStart,chi2,Mode)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else

CommongArg = {'chi2',chi2,'DataEffCor','RunSummary',...
    'ELossFlag','Abdurashitov','exclDataStart',exclDataStart,...
    'SysBudget',0,'DataType','Real','FSDFlag','SAENZ',...
    'fixPar','1 5 6 7 8 9 10 11','NonPoissonScaleFactor',1};

FT1 = MultiRunAnalysis('RunList',RunList1,CommongArg{:});
FT2 = MultiRunAnalysis('RunList',RunList2,CommongArg{:});
FT3 = MultiRunAnalysis('RunList',RunList3,CommongArg{:});
FT4 = MultiRunAnalysis('RunList',RunList4,CommongArg{:});
%% fit
FT1.Fit;
FT2.Fit;
FT3.Fit;
FT4.Fit;

E0 = zeros(4,1); E0Err = zeros(4,1);
E0(1) = FT1.FitResult.par(2); E0Err(1) = FT1.FitResult.err(2);
E0(2) = FT2.FitResult.par(2); E0Err(2) = FT2.FitResult.err(2);
E0(3) = FT3.FitResult.par(2); E0Err(3) = FT3.FitResult.err(2);
E0(4) = FT4.FitResult.par(2); E0Err(4) = FT4.FitResult.err(2);
save(savename,'E0','E0Err','FT1','FT2','FT3','FT4');
end

if contains(Mode,'Scan')
    E0(4) = [];
    E0Err(4) = [];
    x = 1:numel(E0);
else
    x = [FT1.RunData.WGTS_CD_MolPerCm2, FT2.RunData.WGTS_CD_MolPerCm2 ,...
        FT3.RunData.WGTS_CD_MolPerCm2,FT4.RunData.WGTS_CD_MolPerCm2]./FT4.RunData.WGTS_CD_MolPerCm2*100;
end
%% plot

fig1 = figure('Renderer','painters');
set(fig1,'units','centimeters','pos',[0.1, 0.1, 8.4, 4.5]);
plot(linspace(min(x)-20,max(x)+20,numel(x)),zeros(numel(x),1),'-','LineWidth',1,'Color',rgb('Black'));
hold on
e1 = errorbar(x,E0-wmean(E0,1./E0Err.^2),E0Err,'o','MarkerSize',3,'Color',rgb('DarkCyan'),'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
e1.CapSize = 0;
FTpaperFormat;
xlabel(myXlabel);
ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle  (eV)'));

ylim(1.1.*[min(E0-wmean(E0,1./E0Err.^2)-E0Err),max(E0-wmean(E0,1./E0Err.^2)+E0Err)])
set(gca,'FontSize',9);
set(get(gca,'XLabel'),'FontSize',9);
set(get(gca,'YLabel'),'FontSize',9);
switch Mode
    case {'Scan','Scan2'}
        xticks(x);
        xticklabels(myXticks);
        xlim([min(x)-0.2,max(x)+0.2]);
         set(gca,'XMinorTick','off');
    case 'RD'
        xlim([18,105]);
end
%% remove white space
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1)+0.01;
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3)-0.012;
% ax_height = outerpos(4)- ti(2) - ti(4)-0.002;
% ax.Position = [left bottom ax_width ax_height];
set(gca,'LineWidt',1);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(fig1,plotname);







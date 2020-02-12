% Script to check impact of enery loss function on endpoint:
RunList = 'FTpaper';
exclDataStart = 12;
switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7 
        belowE0 = 400;
    case 9
        belowE0 = 200;
    case 12
        belowE0 = 125;
end
ELossAll = {'Aseev','Abdurashitov','KatrinD2'};

savedir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
savename = sprintf('ELossImpactE0_%s_%.0feV.mat',RunList,belowE0);

%if exist([savedir,savename],'file')
%    load([savedir,savename]);
%else
    
    E0Stat = zeros(numel(ELossAll),1);
    E0CM   = zeros(numel(ELossAll),1);
    E0ErrStat = zeros(numel(ELossAll),1);
    E0ErrCM   = zeros(numel(ELossAll),1);
    
    for i=1:numel(ELossAll)
        M = MultiRunAnalysis('RunList',RunList,'exclDataStart',exclDataStart,...
            'ELossFlag',ELossAll{i},'fixPar','1 5 6 7 8 9 10 11');
        M.chi2 = 'chi2Stat';
        M.Fit;
        E0Stat(i)    = M.FitResult.par(2);
        E0ErrStat(i) = M.FitResult.err(2);
        
        M.chi2 = 'chi2CMShape';
        M.ComputeCM;
        M.Fit;
        E0CM(i)    = M.FitResult.par(2);
        E0ErrCM(i) = M.FitResult.err(2);
    end
   
    if ~exist(savedir,'dir')
        system(['mkdir ',savedir])
    end

     save([savedir,savename],'E0Stat','E0ErrStat','E0CM','E0ErrCM');
%end

%% plot
f1 = figure('Name','MultiBarPlot','Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.5]);

x = 1:numel(ELossAll);

eCM = errorbar(x,E0CM+18573.7,E0ErrCM,'o','Color',M.PlotColor,'LineWidth',2,...
    'MarkerFaceColor',M.PlotColor,'MarkerSize',8);
hold on
lCM = plot(linspace(0,max(x)+1,3),mean(E0CM+18573.7).*ones(numel(x),1),'--','LineWidth',2,'Color',M.PlotColor);
eStat =  errorbar(x+0.2,E0Stat+18573.7,E0ErrStat,'o','LineWidth',2,'Color',rgb('GoldenRod'),...
    'MarkerFaceColor',rgb('GoldenRod'),'MarkerSize',8);
lStat = plot(linspace(0,max(x)+1,3),mean(E0Stat+18573.7).*ones(numel(x),1),'--','LineWidth',2,'Color',rgb('GoldenRod'));
xlim([min(x)-0.5,max(x)+0.5]);
xticks([x]);
ytickformat('%.5g')
xticklabels({'Aseev','Abdurashitov','Katrin D2'});
xtickangle=30;
PrettyFigureFormat;
ylabel('E_{0eff.} (eV)');
grid on;
leg = legend([eStat,eCM,lStat,lCM],'stat','stat + sys','mean stat','mean stat + sys');
leg.NumColumns=2;
legend boxoff;
ylim([0.94*min(E0CM-E0ErrCM),1.06*max(E0CM+E0ErrCM)]+18573.7);
% save plot
plotdir =strrep(savedir,'results','plots');
plotname = strrep(savename,'mat','png');

if ~exist(plotdir,'dir')
    system(['mkdir ',plotdir])
end

print(f1,[plotdir,plotname],'-dpng','-r450');






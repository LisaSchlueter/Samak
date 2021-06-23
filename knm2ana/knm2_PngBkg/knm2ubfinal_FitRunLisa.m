% single run fit
% stat + NP only
% nu-mass fix
range   = 40;
freePar = 'E0 Bkg Norm';
chi2    = 'chi2Stat';
SysBudget = 40;
BKG_PtSlope = 3e-06;

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_FitRunList_%.0feV_%s_%s.mat',...
    savedir,range,strrep(freePar,' ',''),chi2);

if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%s.mat',SysBudget));
end

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else   
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2','chi2Stat',...
        'DataType','Real',...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag','StackPixel',...
        'RingMerge','Full',...
        'chi2',chi2,...
        'pullFlag',99,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'DopplerEffectFlag','FSD'...
        'BKG_PtSlope',BKG_PtSlope};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);

    FitResult = A.FitRunList;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A')
end

%% plot
saveplot = 'pdf';
FitResults= A.PlotFitRunList('Parameter','B','YLim',[170 280],'saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','E0','YLim',[-0.7,1],'DisplayStyle','Rel','saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','N','YLim',[0.9,1.1],'saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','pVal','YLim',[-0.2,1.2],'saveplot',saveplot,'HideGaps','OFF');

%%
E0 = A.SingleRun_FitResults.chi2Stat.E0;
B = A.SingleRun_FitResults.chi2Stat.B;
N = A.SingleRun_FitResults.chi2Stat.N;
p = A.SingleRun_FitResults.chi2Stat.pValue;
fprintf('--------------------------------------\n')
fprintf('Run list: %s (%.0f runs) \n',A.RunData.RunName,A.nRuns)
fprintf('<E0> = %.3f eV , std = %.3f eV \n',mean(E0),std(E0));
fprintf('<B>  = %.3f mcps , std = %.3f mcps \n',1e3.*mean(B),1e3.*std(B));
fprintf('<N>  = %.3f       , std = %.3f  \n',1+mean(N),std(N));
fprintf('<pval> = %.2f     , std = %.2f   \n ',mean(p),std(p));
fprintf('--------------------------------------\n')

close all

%% single period
Period = 3;

switch Period
    case 1
        nRuns = 171;
        Start = 1;
        Stop = nRuns;
    case 2
        nRuns = 97;
        Start = 172;
        Stop = Start+nRuns-1;
    case 3
        nRuns = 93;
        Start = 268;
        Stop = 360;
end
fprintf('--------------------------------------\n')
fprintf('Run list Period: %.0f (%.0f runs) \n',Period,nRuns)
fprintf('<E0> = %.3f eV , std = %.3f eV \n',mean(E0(Start:Stop)),std(E0(Start:Stop)));
fprintf('<B>  = %.1f mcps , std = %.1f mcps \n',1e3.*mean(B(Start:Stop)),1e3.*std(B(Start:Stop)));
fprintf('<N>  = %.3f       , std = %.3f  \n',1+mean(N(Start:Stop)),std(N(Start:Stop)));
fprintf('<pval> = %.2f     , std = %.2f   \n ',mean(p(Start:Stop)),std(p(Start:Stop)));
fprintf('--------------------------------------\n')


addPlot = 'ON'; % additional plot
if strcmp(addPlot,'ON')
%% difference

% get old results
E0Old = zeros(A.nRuns,1);
BOld  = zeros(A.nRuns,1);
NOld  = zeros(A.nRuns,1);
pOld  = zeros(A.nRuns,1);

for i=1:1:A.nRuns
    progressbar(i/A.nRuns)
dOld = [getenv('SamakPath'),'tritium-data/fit/Knm2/Uniform_woBPT/','Fit',num2str(A.RunList(i)),'_chi2Stat_39bE0_freeParE0BkgNorm','.mat'];
d = importdata(dOld);
E0Old(i) = d.FitResult.par(2)+A.ModelObj.Q_i;
BOld(i)  = d.FitResult.par(3)+A.ModelObj.BKG_RateSec_i;
NOld(i)  = d.FitResult.par(4);
pOld(i)  = 1-chi2cdf(d.FitResult.chi2min,d.FitResult.dof);
end


%%
plotdir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/plots/FitRunList/'];
MakeDir(plotdir);

fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[0.1, 0.1,1.0,0.5]);

Par = 'N';
switch Par
    case 'E0'
        y = E0-E0Old;
        yStr = sprintf('{\\itE}_0 - {\\itE}_0^{ w/o B_{PT}} (eV)');
    case 'p'
        y = p-pOld;
        yStr = sprintf('{\\itp} - {\\itp}^{ w/o B_{PT}} ');
    case 'B'
         y = (B-BOld).*1e3;
        yStr = sprintf('{\\itB} - {\\itB}^{ w/o B_{PT}} (mcps)');
         case 'N'
         y = (N-NOld);
        yStr = sprintf('{\\itN} - {\\itN}^{ w/o B_{PT}}');
end

LiveTime = hours(A.SingleRunData.StartTimeStamp-A.SingleRunData.StartTimeStamp(1));
plot(LiveTime,y,'o-.','LineWidth',1,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
hold on;
plot(LiveTime,mean(y).*ones(A.nRuns,1),'k-','LineWidth',2);
PrettyFigureFormat
xlabel('Time (hours)');
ylabel(yStr);
leg = legend(sprintf('\\mu = %.3f, \\sigma = %.3f',mean(y),std(y)));
LegAlpha(leg,0.5); leg.FontSize = get(gca,'FontSize');


print(gcf,[plotdir,'FitRunListDiff_',Par,'.png'],'-dpng','-r350');
%% overlay
fig88 = figure('Renderer','painters');
set(fig88,'units','normalized','pos',[0.1, 0.1,1.0,0.5]);
% 
Par = 'B';
switch Par
    case 'E0'
        y1 = E0;
        y2 = E0Old;
        yStr = sprintf('{\\itE}_0 (eV)');
    case 'p'
        y1 = p;
        y2 = pOld;
        yStr = sprintf('{\\itp}-value');
    case 'B'
        y1 = B.*1e3;
        y2 = BOld.*1e3;
        yStr = sprintf('{\\itB} (mcps)');
    case 'N'
        y1 = N+1;
        y2 = NOld+1;
        yStr = sprintf('{\\itN}');
end


pnew = plot(LiveTime,y1,'o-.','LineWidth',1,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
hold on;
pold= plot(LiveTime,y2,'o-.','LineWidth',1,'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
PrettyFigureFormat
xlabel('Time (hours)');
ylabel(yStr);
leg = legend([pnew,pold],sprintf('With {\\itB}_{PT}'),sprintf('Without {\\itB}_{PT}'));
LegAlpha(leg,0.5); leg.FontSize = get(gca,'FontSize');
leg.Location = 'northwest';


print(gcf,[plotdir,'FitRunListOverlay_',Par,'.png'],'-dpng','-r350');
end
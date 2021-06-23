Uniform = 'ON';
MR4 = 'ON';
MR4qU = 'ON';
MR12 = 'ON';
range = 40;
freePar = 'mNu E0 Bkg Norm';
FSDFlag= 'KNM2';
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/BestFit/'];
CommonStrStat =  sprintf('%sknm2ub2_Fit_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2Stat');
CommonStrCM =  sprintf('%sknm2ub2_Fit_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2CMShape');

if strcmp(Uniform,'ON')
    fileStat = [CommonStrStat,sprintf('_StackPixel_%s.mat',FSDFlag)];
    filecm   = [CommonStrCM,sprintf('_StackPixel_%s_SysBudget38.mat',FSDFlag)];
    dUstat    = importdata(fileStat);
    dUcm      = importdata(filecm);
    x = [dUstat.FitResult.par(1),dUcm.FitResult.par(1)];
end

if strcmp(MR4,'ON')
    fileStat = [CommonStrStat,sprintf('_Ring_%s.mat',FSDFlag)];
    filecm   = [CommonStrCM,sprintf('_Ring_%s_SysBudget39.mat',FSDFlag)];
    dMR4stat = importdata(fileStat);
    dMR4cm   = importdata(filecm);
    x = [x,dMR4stat.FitResult.par(1),dMR4cm.FitResult.par(1)];
end

if strcmp(MR4qU,'ON')
    fileStat = strrep([CommonStrStat,sprintf('_Ring_%s.mat',FSDFlag)],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    filecm   = strrep([CommonStrCM,sprintf('_Ring_%s_SysBudget39.mat',FSDFlag)],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    dMR4qUstat = importdata(fileStat);
    dMR4qUcm   = importdata(filecm);
    x = [x,dMR4qUstat.FitResult.par(1),dMR4qUcm.FitResult.par(1)];
end
%%
if strcmp(MR12,'ON')
    fileStat = [CommonStrStat,sprintf('_RingNone_%s.mat',FSDFlag)];
    dMR12stat = importdata(fileStat);
    
    fileStatp  = [CommonStrStat,sprintf('+_RingNone_%s.mat',FSDFlag)];
    dMR12statp = importdata(fileStatp);
    x = [x,dMR12stat.FitResult.par(1),dMR12stat.FitResult.par(1),dMR12statp.FitResult.par(1)];
end

%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.6]);
CommonPlotArg = {'CapSize',0,'LineWidth',2};

pmean = plot(mean(x).*ones(1,10),linspace(0,1e2,10),'LineStyle',':','Color',rgb('Silver'),'LineWidth',2);
hold on;
if strcmp(Uniform,'ON')
    y = [1,1.1];
    eUStat = errorbar(dUstat.FitResult.par(1),y(1),0,0,dUstat.FitResult.errNeg(1),dUstat.FitResult.errPos(1),...
        '.',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',30);
    eUcm = errorbar(dUcm.FitResult.par(1),y(2),0,0,dUcm.FitResult.errNeg(1),dUcm.FitResult.errPos(1),...
        '.',CommonPlotArg{:},'Color',rgb('DodgerBlue'),'MarkerSize',30);
    legStr = {'Uniform'};
else
    y = 0;
    legStr = '';
end

if strcmp(MR4qU,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;
    
    eUStat = errorbar(dMR4qUstat.FitResult.par(1),y(end-1),0,0,dMR4qUstat.FitResult.errNeg(1),dMR4qUstat.FitResult.errPos(1),...
        's',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
    hold on;
    eUcm = errorbar(dMR4qUcm.FitResult.par(1),y(end),0,0,dMR4qUcm.FitResult.errNeg(1),dMR4qUcm.FitResult.errPos(1),...
        's',CommonPlotArg{:},'Color',rgb('DodgerBlue'),'MarkerSize',8,'MarkerFaceColor',rgb('DodgerBlue'));
    legStr = {legStr{:},sprintf('MR-4 (\\DeltaqU)')};
end

if strcmp(MR4,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;    
    eUStat = errorbar(dMR4stat.FitResult.par(1),y(end-1),0,0,dMR4stat.FitResult.errNeg(1),dMR4stat.FitResult.errPos(1),...
        'd',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
    hold on;
    eUcm = errorbar(dMR4cm.FitResult.par(1),y(end),0,0,dMR4cm.FitResult.errNeg(1),dMR4cm.FitResult.errPos(1),...
        'd',CommonPlotArg{:},'Color',rgb('DodgerBlue'),'MarkerSize',8,'MarkerFaceColor',rgb('DodgerBlue'));
    legStr = {legStr{:},'MR-4'};
end


if strcmp(MR12,'ON')
    y(end+1) = y(end)+0.25;% y(end+1) = y(end)+0.1;
    eUStat = errorbar(dMR12stat.FitResult.par(1),y(end),0,0,dMR12stat.FitResult.errNeg(1),dMR12stat.FitResult.errPos(1),...
        '*',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
    eUStatp = errorbar(dMR12statp.FitResult.par(1),y(end)+0.1,0,0,dMR12statp.FitResult.errNeg(1),dMR12statp.FitResult.errPos(1),...
        '*',CommonPlotArg{:},'Color',rgb('SkyBlue'),'MarkerSize',8,'MarkerFaceColor',rgb('DodgerBlue'));
    legStr = {legStr{:},'MR-12'};
end

% tweak plot for legend
pcm   = plot(0,1e2,'LineWidth',2,'Color',rgb('DodgerBlue'));
pstat =  plot(0,1e2,'LineWidth',2,'Color',rgb('Orange'));
%
PrettyFigureFormat('FontSize',22)
ylim([y(1)-0.1,y(end)+0.35]);%+0.18]);
yticks([y(1)+0.05,y(3)+0.05,y(5)+0.05,y(end)+0.05]); set(gca,'YMinorTick','off');
yticklabels(legStr);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
leg = legend([pstat,pcm],'Stat. only','Stat. and syst.','EdgeColor',rgb('Silver'),'Location','northwest');
%% save
plotdir = strrep(savedir,'results/BestFit','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm2ub2_FitResultOverview_mNuSq_%s.png',plotdir,FSDFlag);
print(plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname)
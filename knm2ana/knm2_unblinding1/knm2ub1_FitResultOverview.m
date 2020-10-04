Uniform = 'ON';
MR4 = 'ON';
MR12 = 'ON';
range = 40;
freePar = 'mNu E0 Bkg Norm';
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/BestFit/'];
CommonStrStat =  sprintf('%sknm2ub1_Fit_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2Stat');
CommonStrCM =  sprintf('%sknm2ub1_Fit_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2CMShape');

if strcmp(Uniform,'ON')
    fileStat = [CommonStrStat,'_StackPixel.mat'];
    filecm   = [CommonStrCM,'_StackPixel_SysBudget38.mat'];
    dUstat    = importdata(fileStat);
    dUcm      = importdata(filecm);
    x = [dUstat.FitResult.par(1),dUcm.FitResult.par(1)];
end

if strcmp(MR4,'ON')
    fileStat = [CommonStrStat,'_Ring.mat'];
    filecm   = [CommonStrCM,'_Ring_SysBudget39.mat'];
    dMR4stat = importdata(fileStat);
    dMR4cm   = importdata(filecm);
    x = [x,dMR4stat.FitResult.par(1),dMR4cm.FitResult.par(1)];
end
%%
if strcmp(MR12,'ON')
    fileStat = [CommonStrStat,'_RingNone.mat'];
    dMR12stat = importdata(fileStat);
    x = [x,dMR12stat.FitResult.par(1),dMR12stat.FitResult.par(1)];
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
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;
    eUStat = errorbar(dMR12stat.FitResult.par(1),y(end-1),0,0,dMR12stat.FitResult.errNeg(1),dMR12stat.FitResult.errPos(1),...
        'd',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
    legStr = {legStr{:},'MR-12'};
end

% tweak plot for legend
pcm   = plot(0,1e2,'LineWidth',2,'Color',rgb('DodgerBlue'));
pstat =  plot(0,1e2,'LineWidth',2,'Color',rgb('Orange'));
%%
PrettyFigureFormat('FontSize',22)
ylim([y(1)-0.1,y(end)+0.1]);%+0.18]);
yticks([y(1)+0.1,y(3)+0.1,y(5)]); set(gca,'YMinorTick','off');
yticklabels(legStr);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
leg = legend([pstat,pcm],'Stat. only','Stat. and syst.','EdgeColor',rgb('Silver'),'Location','north');
 
plotdir = strrep(savedir,'results/BestFit','plots');
plotname = sprintf('%sknm2ub1_FitResultOverview_mNuSq.png',plotdir);
print(plotname,'-dpng','-r350');
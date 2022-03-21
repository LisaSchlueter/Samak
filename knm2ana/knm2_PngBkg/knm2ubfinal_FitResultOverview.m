% comparison plot of unblinded fit results and upper limits

Uniform = 'ON';
MR4 = 'OFF';
MR4qU = 'ON';
MR12qU = 'ON';
range = 40;
freePar = 'mNu E0 Bkg Norm';
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
CommonStrStat =  sprintf('%sknm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2Stat');
CommonStrCM =  sprintf('%sknm2ubfinal_Fit_Bpng-3.0mucpsPers_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2CMShape');

if strcmp(Uniform,'ON')
    fileStat = [CommonStrStat,'_StackPixel_KNM2.mat'];
    filecm   = [CommonStrCM,'_StackPixel_KNM2_SysBudget40.mat'];
    dUstat    = importdata(fileStat);
    dUcm      = importdata(filecm);
    x      = [dUstat.FitResult.par(1),dUcm.FitResult.par(1)];
end

if strcmp(MR4,'ON')
    fileStat = [CommonStrStat,'_Ring_KNM2.mat'];
    filecm   = [CommonStrCM,'_Ring_KNM2_SysBudget39.mat'];
    dMR4stat = importdata(fileStat);
    dMR4cm   = importdata(filecm);
    x = [x,dMR4stat.FitResult.par(1),dMR4cm.FitResult.par(1)];
end

if strcmp(MR4qU,'ON')
    fileStat = strrep([CommonStrStat,'_Ring_KNM2.mat'],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    filecm   = strrep([CommonStrCM,'_Ring_KNM2_SysBudget41.mat'],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    dMR4qUstat = importdata(fileStat);
    dMR4qUcm   = importdata(filecm);
    x = [x,dMR4qUstat.FitResult.par(1),dMR4qUcm.FitResult.par(1)];
end

if strcmp(MR12qU,'ON')
     fileStat = strrep([CommonStrStat,'_RingNone_KNM2_pull6.mat'],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    dMR12stat = importdata(fileStat);
     fileCM = strrep([CommonStrCM,'_RingNone_KNM2_SysBudget40_pull6.mat'],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    dMR12cm = importdata(fileCM);
    
    x = [x,dMR12stat.FitResult.par(1),dMR12cm.FitResult.par(1)];
end

%% load upper limit as a function of best fit
savenameFC = sprintf('%sknm2ub2_mNuLimits_%s.mat',savedir,'FC');
savenameLT = sprintf('%sknm2ub2_mNuLimits_%s.mat',savedir,'LT');
dFC = importdata(savenameFC);
dLT = importdata(savenameLT);

UpperLimFC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,x,'spline'));
UpperLimLT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,x,'spline'));
%
%% plot
LocalFontSize = 18;
UpperLim = 'ON';
close all
if strcmp(UpperLim,'ON')
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.45]);
    s1 = subplot(1,3,1:2);
else
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.6]);
end

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
    legStr = {legStr{:},sprintf('MR-4')};
end

if strcmp(MR4,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;    
    eUStat = errorbar(dMR4stat.FitResult.par(1),y(end-1),0,0,dMR4stat.FitResult.errNeg(1),dMR4stat.FitResult.errPos(1),...
        'd',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
    hold on;
    eUcm = errorbar(dMR4cm.FitResult.par(1),y(end),0,0,dMR4cm.FitResult.errNeg(1),dMR4cm.FitResult.errPos(1),...
        'd',CommonPlotArg{:},'Color',rgb('DodgerBlue'),'MarkerSize',8,'MarkerFaceColor',rgb('DodgerBlue'));
    legStr = {legStr{:},'MR-4 (no qU)'};
end


if strcmp(MR12qU,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;
    eUStat = errorbar(dMR12stat.FitResult.par(1),y(end-1),0,0,dMR12stat.FitResult.errNeg(1),dMR12stat.FitResult.errPos(1),...
        '*',CommonPlotArg{:},'Color',rgb('Orange'),'MarkerSize',8,'MarkerFaceColor',rgb('Orange'));
     eUcm = errorbar(dMR12cm.FitResult.par(1),y(end),0,0,dMR12cm.FitResult.errNeg(1),dMR12cm.FitResult.errPos(1),...
        '*',CommonPlotArg{:},'Color',rgb('DodgerBlue'),'MarkerSize',8,'MarkerFaceColor',rgb('DodgerBlue'));
    legStr = {legStr{:},'MR-12'};
end

% tweak plot for legend
pcm   = plot(0,1e2,'LineWidth',2,'Color',rgb('DodgerBlue'));
pstat =  plot(0,1e2,'LineWidth',2,'Color',rgb('Orange'));
%

ylim([y(1)-0.1,y(end)+0.22]);%+0.18]);
if strcmp(MR12qU,'ON')
    if strcmp(MR4,'ON')
    yticks([y(1)+0.05,y(3)+0.05,y(5)+0.05]);
    else
      yticks([y(1)+0.05,y(3)+0.05,y(end)]);  
    end
else
    yticks([y(1)+0.05,y(3)+0.05]);
end
set(gca,'YMinorTick','off');
yticklabels(legStr);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^{ 2})'));
leg = legend([pstat,pcm],'Stat. only','Stat. and syst.','EdgeColor',rgb('Silver'),'Location','northwest');
PrettyLegendFormat(leg);
PrettyFigureFormat('FontSize',LocalFontSize)
xlim([-0.1,0.65])
ax1 = gca;
ax1.YAxis.FontSize = ax1.XLabel.FontSize;
leg.FontSize = get(gca,'FontSize')+2;

if strcmp(UpperLim,'ON')
    s2 = subplot(1,3,3);
    PArg = {'MarkerSize',8,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'),'LineWidth',2};
    plot(UpperLimFC(2),y(2),'o',PArg{:});
    hold on;
    plot(UpperLimFC(4),y(4),'s',PArg{:});
    plot(UpperLimFC(6),y(6),'*',PArg{:});
    
    % UpLT = plot(UpperLimLT(1:2:end),y(1:2:end),'d','MarkerSize',6,'Color',rgb('Orange'));
    % plot(UpperLimLT(2:2:end),y(2:2:end),'x','MarkerSize',6,'Color',rgb('DodgerBlue'));
    linkaxes([s1,s2],'y');
    
    xlabel(sprintf('Upper Limit in {\\itm}_\\nu (eV)\n at 90%% C.L.'));
     PrettyFigureFormat('FontSize',LocalFontSize)
    ylim([y(1)-0.1,y(end)+0.22]);%+0.18]);
    if strcmp(MR12qU,'ON')
        if strcmp(MR4,'ON')
            yticks([y(1)+0.05,y(3)+0.05,y(5)+0.05]);
        else
            yticks([y(1)+0.05,y(3)+0.05,y(end)]);
        end
    else
        yticks([y(1)+0.05,y(3)+0.05]);
    end
    yticklabels([])
    set(gca,'YMinorTick','off'); 
   
    xlim([0.886 0.915])
    xticks([0.89:0.01:0.91])
   
end

ax2 = gca;
ax1.Position(2) = 0.18;
ax2.Position(2) = 0.18;
ax2.Position(1) = 0.635;
ax2.Position(3) = 0.2;
ax1.XLabel.Position(2) = ax2.XLabel.Position(2);

%
plotdir = strrep(savedir,'results','plots');
plotname = sprintf('%sknm2ubfinal_FitResultOverview_mNuSq_UpLim%s.pdf',plotdir,UpperLim);
export_fig(plotname);
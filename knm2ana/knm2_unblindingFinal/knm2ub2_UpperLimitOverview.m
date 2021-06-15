%% load upper limits

savedir  = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/'];
savenameFC = sprintf('%sknm2ub2_mNuLimits_%s.mat',savedir,'FC');
savenameLT = sprintf('%sknm2ub2_mNuLimits_%s.mat',savedir,'LT');
dFC = importdata(savenameFC);
dLT = importdata(savenameLT);

%mNuLim_FC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,x,'spline'));
%mNuLim_LT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,x,'spline'));
%%

Uniform = 'ON';
MR4 = 'ON';
MR4qU = 'ON';
MR12 = 'ON';
range = 40;
freePar = 'mNu E0 Bkg Norm';
FSDFlag= 'KNM2';
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/BestFit/'];
CommonStrStat =  sprintf('%sknm2ub2_Fit_Real_%.0feV_%s_%s',savedir,range,strrep(freePar,' ',''),'chi2Stat');

if strcmp(Uniform,'ON')
    filecm   = [CommonStrCM,sprintf('_StackPixel_%s_SysBudget38.mat',FSDFlag)];
    dUcm      = importdata(filecm);
    dUcm.UpperLimFC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,dUcm.FitResult.par(1),'spline'));
    dUcm.UpperLimLT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,dUcm.FitResult.par(1),'spline'));
end

if strcmp(MR4,'ON')
    filecm   = [CommonStrCM,sprintf('_Ring_%s_SysBudget39.mat',FSDFlag)];
    dMR4cm   = importdata(filecm);
    dMR4cm.UpperLimFC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,dMR4cm.FitResult.par(1),'spline'));
    dMR4cm.UpperLimLT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,dMR4cm.FitResult.par(1),'spline'));
end

if strcmp(MR4qU,'ON')
    filecm   = strrep([CommonStrCM,sprintf('_Ring_%s_SysBudget39.mat',FSDFlag)],strrep(freePar,' ',''),[strrep(freePar,' ',''),'qU']);
    dMR4qUcm   = importdata(filecm);
    dMR4qUcm.UpperLimFC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,dMR4qUcm.FitResult.par(1),'spline'));
    dMR4qUcm.UpperLimLT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,dMR4qUcm.FitResult.par(1),'spline'));
end
%
if strcmp(MR12,'ON')
    fileStat = [CommonStrStat,sprintf('_RingNone_%s.mat',FSDFlag)];
    dMR12stat = importdata(fileStat);
    dMR12stat.UpperLimFC = sqrt(interp1(dFC.mNuMeasured_v,dFC.mNuSqLimit_v,dMR12stat.FitResult.par(1),'spline'));
    dMR12stat.UpperLimLT = sqrt(interp1(dLT.mNuMeasured_v,dLT.mNuSqLimit_v,dMR12stat.FitResult.par(1),'spline'));
end



%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.6]);
CommonPlotArg = {'LineWidth',2};

%pmean = plot(mean(x).*ones(1,10),linspace(0,1e2,10),'LineStyle',':','Color',rgb('Silver'),'LineWidth',2);
if strcmp(Uniform,'ON')
    y = [1,1.1];
    plot(dUcm.UpperLimFC,y(1),'.',CommonPlotArg{:},'Color',rgb('Tomato'),'MarkerSize',30);
    hold on;
    plot(dUcm.UpperLimLT,y(2), '.',CommonPlotArg{:},'Color',rgb('Turquoise'),'MarkerSize',30);
    legStr = {'Uniform'};
else
    y = 0;
    legStr = '';
end

if strcmp(MR4qU,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;
    
    plot(dMR4qUcm.UpperLimFC,y(end-1),...
        's',CommonPlotArg{:},'Color',rgb('Tomato'),'MarkerSize',8,'MarkerFaceColor',rgb('Tomato'));
    hold on;
   plot(dMR4qUcm.UpperLimLT,y(end),...
        's',CommonPlotArg{:},'Color',rgb('Turquoise'),'MarkerSize',8,'MarkerFaceColor',rgb('Turquoise'));
    legStr = {legStr{:},sprintf('MR-4 (\\DeltaqU)')};
end

if strcmp(MR4,'ON')
    y(end+1) = y(end)+0.25; y(end+1) = y(end)+0.1;    
    plot(dMR4cm.UpperLimFC,y(end-1),...
        'd',CommonPlotArg{:},'Color',rgb('Tomato'),'MarkerSize',8,'MarkerFaceColor',rgb('Tomato'));
    hold on;
    plot(dMR4cm.UpperLimLT,y(end),...
        'd',CommonPlotArg{:},'Color',rgb('Turquoise'),'MarkerSize',8,'MarkerFaceColor',rgb('Turquoise'));
    legStr = {legStr{:},'MR-4'};
end


if strcmp(MR12,'ON')
    y(end+1) = y(end)+0.25;% y(end+1) = y(end)+0.1;
    plot(dMR12stat.UpperLimFC,y(end),...
        '*',CommonPlotArg{:},'Color',rgb('Tomato'),'MarkerSize',8,'MarkerFaceColor',rgb('Tomato'));
    plot(dMR12stat.UpperLimLT,y(end)+0.1,...
        '*',CommonPlotArg{:},'Color',rgb('Turquoise'),'MarkerSize',8,'MarkerFaceColor',rgb('Turquoise'));
    legStr = {legStr{:},'MR-12'};
end

% tweak plot for legend
pLT  = plot(0,1e2,'LineWidth',2,'Color',rgb('Turquoise'));
pFC =  plot(0,1e2,'LineWidth',2,'Color',rgb('Tomato'));
%
PrettyFigureFormat('FontSize',22)
ylim([y(1)-0.1,y(end)+0.22]);%+0.18]);
yticks([y(1)+0.05,y(3)+0.05,y(5)+0.05,y(end)]); set(gca,'YMinorTick','off');
yticklabels(legStr);
xlim([0.75 0.9])
grid on;
xlabel(sprintf('Upper limit on {\\itm}_\\nu at 90%% C.L. (eV)'));
leg = legend([pLT,pFC],'Lokhov-Tkachov','Feldman-Cousins','EdgeColor',rgb('Silver'),'Location','northwest');
 %% save
plotdir = strrep(savedir,'results/BestFit','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm2ub2_UpperLimOverview_mNuSq_%s.png',plotdir,FSDFlag);
print(plotname,'-dpng','-r350');
fprintf('save plot to %s \n',plotname)
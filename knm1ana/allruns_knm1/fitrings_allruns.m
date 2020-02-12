range = 90;
if range==90
elseif range==40
end
%% -175mV rear wall
M175 = MultiRunAnalysis('RunList','KNM1_175mvRW');
R175 = RingAnalysis('RunAnaObj',M175,'RingList',1:12);
R175.FitRings('SaveResult','ON','RecomputeFlag','OFF');
R175.PlotFits('SavePlot','ON');

%% -183mV rear wall
M183 = MultiRunAnalysis('RunList','KNM1_m183mvRW');
R183 = RingAnalysis('RunAnaObj',M183,'RingList',1:12);
R183.FitRings('SaveResult','ON','RecomputeFlag','OFF');
R183.PlotFits('SavePlot','ON');

%% -149mV rear wall
M149 = MultiRunAnalysis('RunList','KNM1','RingMerge','Full','exclDataStart',13,'DataType','Twin');
R149 = RingAnalysis('RunAnaObj',M149,'RingList',1:4);
R149.FitRings('SaveResult','OFF','RecomputeFlag','ON');
R149.PlotFits('SavePlot','OFF','PlotPar',2,'Blind','OFF');

%% do linear fits all together
PlotPar = 2;
meanPar175 =wmean(R175.FitResult.par(:,PlotPar),1./R175.FitResult.err(:,PlotPar).^2);
meanPar183 =wmean(R183.FitResult.par(:,PlotPar),1./R183.FitResult.err(:,PlotPar).^2);
meanPar149 =wmean(R149.FitResult.par(:,PlotPar),1./R149.FitResult.err(:,PlotPar).^2);
[linFitpar175, linFiterr175, linFitchi2min175,linFitdof175] =...
    linFit(R175.RingList',R175.FitResult.par(:,PlotPar)-meanPar175,R175.FitResult.err(:,PlotPar));
[linFitpar183, linFiterr183, linFitchi2min183,linFitdof183] =...
    linFit(R183.RingList',R183.FitResult.par(:,PlotPar)-meanPar183,R183.FitResult.err(:,PlotPar));
[linFitpar149, linFiterr149, linFitchi2min149,linFitdof149] =...
    linFit(R149.RingList',R149.FitResult.par(:,PlotPar)-meanPar149,R149.FitResult.err(:,PlotPar));
%% plot
plotOpt = 'all';
linFitFlag = 'ON';
fig2 = figure('Renderer','opengl');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);

c175 = rgb('Orange'); c183 = rgb('CadetBlue'); c149 = rgb('FireBrick');

e175 = errorbar(R175.RingList,R175.FitResult.par(:,PlotPar)-meanPar175,R175.FitResult.err(:,PlotPar),'o',...
    'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',c175,'MarkerFaceColor',c175);
hold on;
plot(linspace(0.5,12.5,R175.nRings),zeros(R175.nRings,1),'--','Color',rgb('Black'),'LineWidth',2);
l175 = plot(R175.RingList,linFitpar175(1).*R175.RingList+linFitpar175(2),'-','Color',c175,'LineWidth',3);
linFitleg175 =  sprintf('linear fit slope: (%.1f \\pm %.1f) meV @ \\chi2 = %.1f / %.0f dof',linFitpar175(1)*1e3,linFiterr175(1)*1e3,linFitchi2min175,linFitdof175);

e183 = errorbar(R183.RingList,R183.FitResult.par(:,PlotPar)-meanPar183,R183.FitResult.err(:,PlotPar),'o',...
    'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',c183,'MarkerFaceColor',c183);
l183 = plot(R183.RingList,linFitpar183(1).*R183.RingList+linFitpar183(2),'-','Color',c183,'LineWidth',3);
linFitleg183 =  sprintf('linear fit slope: (%.1f \\pm %.1f) meV @ \\chi2 = %.1f / %.0f dof',linFitpar183(1)*1e3,linFiterr183(1)*1e3,linFitchi2min183,linFitdof183);


e149 = errorbar(R149.RingList,R149.FitResult.par(:,PlotPar)-meanPar149,R149.FitResult.err(:,PlotPar),'o',...
    'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',c149,'MarkerFaceColor',c149);
l149 = plot(R149.RingList,linFitpar149(1).*R149.RingList+linFitpar149(2),'-','Color',c149,'LineWidth',3);
linFitleg149 =  sprintf('linear fit slope: (%.1f \\pm %.1f) meV @ \\chi2 = %.1f / %.0f dof',linFitpar149(1)*1e3,linFiterr149(1)*1e3,linFitchi2min149,linFitdof149);

leg = legend([e175,e183,e149,l175,l183,l149],...
    'RW  175mV','RW -183mV','RW -149mV',linFitleg175,linFitleg183,linFitleg149);
leg.NumColumns = 2; legend boxoff
leg.Location = 'northwest';
% style
hold off;
PrettyFigureFormat;
if PlotPar==1
    ylabel('m^2_\nu - <m^2_\nu> (eV^2)');
elseif PlotPar==2
    ylabel('E_0 - <E_0> (eV)');
elseif PlotPar==3
    ylabel('B - <B> (cps)');
elseif  PlotPar==4
    ylabel('N - <N> ');
end
xlabel('ring')
set(gca,'FontSize',20);

xlim([0.5,12.5])


if PlotPar==2
    Parlabel='E0';
elseif PlotPar==3
    Parlabel='B';
elseif PlotPar==4
    Parlabel='N';
end
savename = sprintf('./plots/%sringwise_all_%.0frange.png',Parlabel,range);
print(savename,'-dpng','-r500');
%%
%M149.Fit;M183.Fit; M175.Fit;
fig1 = figure('Renderer','opengl');
set(fig1,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);
errorbar([-183,-149,175],1e3.*([M183.FitResult.par(2),M149.FitResult.par(2),M175.FitResult.par(2)]-M149.FitResult.par(2)),1e3.*[M183.FitResult.err(2),M149.FitResult.err(2),M175.FitResult.err(2)],...
    '--o','LineWidth',3,'Color',rgb('IndianRed'),'MarkerSize',10);
xlabel('rear wall voltage (mV)');
ylabel('rel. endpoint (meV)');
grid on;
PrettyFigureFormat;
set(gca,'FontSize',18)

savename = sprintf('./plots/%sE0shift_all_%.0frange.png',Parlabel,range);
print(savename,'-dpng','-r500');




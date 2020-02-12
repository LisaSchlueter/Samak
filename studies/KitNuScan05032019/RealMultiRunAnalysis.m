%
M=MultiRunAnalysis('DataType','Real','RunList','KITNUSCAN','exclDataStart',2,'fixPar',' 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
M.Fit;

%
Mu=MultiRunAnalysis('DataType','Real','RunList','KITNUSCANup','exclDataStart',2,'fixPar','1 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
Mu.Fit;
Mu.PlotFit('saveplot','png');
qU_up = Mu.RunData.qU(30:end);
B_up  = Mu.RunData.TBDIS(30:end)./(Mu.RunData.qUfrac(30:end).*Mu.RunData.TimeSec);
Be_up  = Mu.RunData.TBDISE(30:end)./(Mu.RunData.qUfrac(30:end).*Mu.RunData.TimeSec);

%
Md=MultiRunAnalysis('DataType','Real','RunList','KITNUSCANdown','exclDataStart',2,'fixPar','1 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
Md.Fit;
Md.PlotFit('saveplot','png');
qU_down = Md.RunData.qU(30:end);
B_down  = Md.RunData.TBDIS(30:end)./(Md.RunData.qUfrac(30:end).*Md.RunData.TimeSec);
Be_down  = Md.RunData.TBDISE(30:end)./(Md.RunData.qUfrac(30:end).*Md.RunData.TimeSec);

%% Plot Data
figure(070320191)
% show the baseline
plot([min(qU_up) max(qU_up)], [M.ModelObj.BKG_RateSec M.ModelObj.BKG_RateSec], 'k-', 'linewidth', 1);
hold on
errorbar(qU_up,B_up,Be_up,'markersize', 12, 'Marker', 'o', 'markerfacecolor', rgb('DarkBlue'),'markeredgecolor',rgb('DarkBlue'),'LineWidth',2);
errorbar(qU_down,B_down,Be_down,'markersize', 12, 'Marker', '^', 'markerfacecolor', rgb('Amethyst'),'markeredgecolor',rgb('Amethyst'),'LineWidth',2);
hold off
PrettyFigureFormat
xlabel('qU (eV)');
ylabel('cps');

%% Kolmogorov-Smirnov Test
[h,p,k] = kstest2(B_up,B_down);
figure(070320192)
h1=cdfplot(B_up);h1.LineWidth=2;
hold on
h2=cdfplot(B_down);h2.LineWidth=2;
legend([h1, h2],'Up Scans','Down Scans','Location','best');
PrettyFigureFormat
hold off
xlabel('background (cps)');
ylabel('cumulative pdf');
if h==0
    title('Up/Down Background Distribution Compatible (95%CL)');
else
    title('Up/Down Background Distribution Different (95%CL)');
end

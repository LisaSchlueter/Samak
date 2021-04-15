% Test of Wilk's theorem (coverage)
%  DeltaChi2 distribution (best fit - null chi2)
% cumulative pdf with critical delta chi2 for 95%CL.

NrandMC = 606;
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file')
    load(savefile);
    fprintf('load file from %s \n',savefile);
else
     fprintf('file does not exist: %s \n',savefile);
     return
end

%%

%% DeltaChi2 distribution (best fit - null chi2)

PlotDeltaChi2 = sort(chi2_delta);
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);

% calculate 95 quantile: interpolation
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2CrApprox = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'lin');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));
x = linspace(0,max(PlotDeltaChi2),1e3);
DeltaChi2CDFTheo = chi2cdf(x,2);
xInter = linspace(0,10,1e2);
DeltaChi2CrTheo = interp1(chi2cdf(xInter,2),xInter,0.95,'spline');

%% plot cdf
GetFigure;
p95 = plot(linspace(0,20,10),0.95*ones(10,1),'-','LineWidth',1.5,'Color',rgb('Silver'));
hold on;
pchi2Theo = plot(x,DeltaChi2CDFTheo,'-','LineWidth',2.5,'Color',rgb('Orange'));
pchi2 = plot(PlotDeltaChi2,DeltaChi2CDF,'-.','LineWidth',2.5,'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('\\Delta \\chi^2'));
ylabel(sprintf('Cumulative probability'));
%xlim([0 10]);
ylim([-0.02 1.02]);
legend([p95,pchi2Theo,pchi2,],sprintf('95%% quantile'),...
    sprintf('\\chi^2 distribution for 2 dof          \\Delta\\chi^2_{crit.} = %.2f',DeltaChi2CrTheo),...
    sprintf('Empirical cdf (%.0f samples) \\Delta\\chi^2_{crit.} = %.2f',numel(PlotDeltaChi2),DeltaChi2CrApprox),...
    'EdgeColor',rgb('Silver'),'Location','southeast');
t = title(sprintf('\\Delta\\chi^2 = \\chi^2_{null} - \\chi^2_{min} '),'FontWeight','normal','FontSize',get(gca,'FontSize'));
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');

%% save
% plotname = strrep(strrep(savefile,'results','plots'),'.mat','_DeltaChi2Crit.png');
% print(gcf,plotname,'-dpng','-r450');
% fprintf('save plot to %s \n',plotname);
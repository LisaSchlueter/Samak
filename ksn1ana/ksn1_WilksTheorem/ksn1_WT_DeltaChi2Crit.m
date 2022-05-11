range = 40;%
SavePlt = 'ON';
chi2Str = 'chi2CMShape';
InterpMode = 'lin';
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_WilksTheorem/results/'];
MakeDir(savedir);
savename = sprintf('%sksn1_WilksTheorem_%.0frange_%s_%s.mat',savedir,range,chi2Str,InterpMode);

if exist(savename,'file')
    d = importdata(savename);
    fprintf('load %s\n',savename);
else
    fprintf(2,'file not found, run ksn1_WT_MergeFiles.m \n');
    return
end

chi2_delta = d.chi2min_null-d.chi2min_bf;
%% DeltaChi2 distribution (best fit - null chi2)
dof = 2;
chi2_delta(chi2_delta<0) = 0;
PlotDeltaChi2 = sort(chi2_delta);%unique(chi2_delta);%
DeltaChi2CDF = arrayfun(@(x) sum(PlotDeltaChi2<=x)./numel(PlotDeltaChi2),PlotDeltaChi2);

% calculate 95 quantile: interpolation
[DeltaChi2CDFquantile,ia] = unique(DeltaChi2CDF);
DeltaChi2CrApprox = interp1(DeltaChi2CDFquantile,PlotDeltaChi2(ia),0.95,'lin');%quantile(PlotDeltaChi2,0.95);% PlotDeltaChi2(find(abs(DeltaChi2CDF-0.95)==min(abs(DeltaChi2CDF-0.95)),1));
x = linspace(0,max(PlotDeltaChi2),1e3);
DeltaChi2CDFTheo = chi2cdf(x,dof);
xInter = linspace(0,10,1e2);
DeltaChi2CrTheo = interp1(chi2cdf(xInter,dof),xInter,0.95,'spline');


%% plot cdf
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
p95 = plot(linspace(0,20,10),0.95*ones(10,1),'-','LineWidth',3,'Color',rgb('Silver'));
hold on;
pchi2Theo = plot(x,DeltaChi2CDFTheo,'-','LineWidth',3,'Color',rgb('Orange'));
pchi2 = plot(PlotDeltaChi2,DeltaChi2CDF,'-.','LineWidth',3,'Color',rgb('DodgerBlue'));

xlabel(sprintf('\\Delta \\chi^2'));
ylabel(sprintf('Cumulative probability'));
PrettyFigureFormat('FontSize',22);
xlim([0 10]);
ylim([0 1.02]);

leg = legend([p95,pchi2Theo,pchi2,],sprintf('95%% quantile'),...
    sprintf('Chi-squared distribution for 2 dof:     \\Delta\\chi^2_{crit.} = %.2f',DeltaChi2CrTheo),...
    sprintf('Empirical KNM1 cdf (%.0f samples): \\Delta\\chi^2_{crit.} = %.2f',numel(PlotDeltaChi2),DeltaChi2CrApprox),...
    'EdgeColor',rgb('Silver'),'Location','southeast');

% 
% leg = legend([p95,pchi2Theo,pchi2,],sprintf('95%% quantile'),...
%     sprintf('\\chi^2 distribution for 2 dof          \\Delta\\chi^2_{crit.} = %.2f',DeltaChi2CrTheo),...
%     sprintf('Empirical cdf (%.0f samples) \\Delta\\chi^2_{crit.} = %.2f',numel(PlotDeltaChi2),DeltaChi2CrApprox),...
%     'EdgeColor',rgb('Silver'),'Location','southeast');
leg.FontSize = get(gca,'FontSize')+2;
PrettyLegendFormat(leg);
%t = title(sprintf('\\Delta\\chi^2 = \\chi^2_{null} - \\chi^2_{min} '),'FontWeight','normal','FontSize',get(gca,'FontSize'));
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');

%% save
%  plotname = strrep(strrep(savefile,'results','plots'),'.mat','_DeltaChi2Crit.png');
%  print(gcf,plotname,'-dpng','-r450');
%  fprintf('save plot to %s \n',plotname);
 
 %% KS Test
CDFTheo = chi2cdf(PlotDeltaChi2,2); 
 [h,p,ksstat,cv] = kstest(PlotDeltaChi2,'CDF',[PlotDeltaChi2,CDFTheo]);
 %%
 pltdir= strrep(savedir,'results','plots');
 pltname = sprintf('%sksn1_WT_DeltaChi2CritH0.pdf',pltdir);
   export_fig(pltname);
fprintf('save plot to %s \n',pltname);
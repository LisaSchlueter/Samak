% stupid plot script to plot neyman upper limit

%get stat only values
M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11');
S = RunSensitivity('RunAnaObj',M);
[UpperLimitStat,mNuSq_Quantil90Stat,mNuSq_tStat]  = ...
    S.ComputeUpperLimit('nSamples',1000,'mNuSq_t',[0,0.1,0.2,0.3,0.7,2,5:0.5:6]);
mNuMeasuredStat = min(mNuSq_Quantil90Stat):0.05:3;
MyUpperLimitsStat = interp1(mNuSq_Quantil90Stat,mNuSq_tStat,mNuMeasuredStat);

%% get stat + sys values
M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11',...
    'SysBudget',2,'chi2','chi2CMShape');
S = RunSensitivity('RunAnaObj',M);


[UpperLimit,mNuSq_Quantil90,mNuSq_t]  = ...
    S.ComputeUpperLimit('nSamples',1000,'mNuSq_t',[0,0.5,1,1.2,2,3,4]);%1.2:-0.2:0.2);%[0,0.5,1,1.5,2,3,0.2:0.2:1.2]);
mNuMeasuredCM= min(mNuSq_Quantil90):0.05:3;
MyUpperLimitsCM = interp1(mNuSq_Quantil90,mNuSq_t,mNuMeasuredCM);
%% plot both Upper Limit
f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

pstat  = plot(mNuMeasuredStat,sqrt(MyUpperLimitsStat),'LineWidth',4,'Color',rgb('IndianRed'));
hold on;
psys  = plot(mNuMeasuredCM,sqrt(MyUpperLimitsCM),'LineWidth',4,'Color',rgb('SteelBlue'));
PrettyFigureFormat;
xlabel(sprintf('m^2_\\nu measured (eV^2)'));
ylabel(sprintf('m_\\nu  upper limit (eV) %.0f%% C.L.',S.ConfLevel*100));
leg = legend('stat','stat + sys'); legend boxoff;
leg.Location='northwest';
leg.Title.String = sprintf('%.0f eV range',S.GetRange);
grid on;
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/plots/'];
leg.Location = 'northwest';
savename = [savedir,'Neyman_CL_stat_statsys.png'];
set(gca,'FontSize',22);
xlim([-1.6,2.2]);
print(f1,savename,'-dpng','-r450')
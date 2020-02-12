if ~exist('M','var')
    M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11',...
        'SysBudget',22,'chi2','chi2CMShape');
end

if ~exist('S','var')
    S = RunSensitivity('RunAnaObj',M);
end

[UpperLimit,mNuSq_Quantil90,mNuSq_t]  = ...
    S.ComputeUpperLimit('nSamples',1000,'mNuSq_t',[0.95^2,0]);%[0,0.5,1.2]1.2:-0.2:0.2);%[0,0.5,1,1.5,2,3,0.2:0.2:1.2]);

%%%
% plot Upper Limit
f1 = figure('Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

mNuMeasured = min(mNuSq_Quantil90):0.05:3;
MyUpperLimits = interp1(mNuSq_Quantil90,mNuSq_t,mNuMeasured);
p1  = plot(mNuMeasured,sqrt(MyUpperLimits),'LineWidth',4,'Color',rgb('IndianRed'));
PrettyFigureFormat;
xlabel(sprintf('m^2_\\nu measured (eV^2)'));
ylabel(sprintf('m_\\nu  upper limit (eV) %.0f%% C.L.',S.ConfLevel*100));
leg = legend(sprintf('stat + sys \n%.0f eV range',S.GetRange)); legend boxoff;
leg.Location='northwest';
grid on;
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/plots/'];
leg.Location = 'northwest';
savename = [savedir,'Neyman_statSysCL.png'];
set(gca,'FontSize',22);
print(f1,savename,'-dpng','-r450')
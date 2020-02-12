
M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','fixPar','5 6 7 8 9 10 11 12',...
    'exclDataStart',14,'FSD','BlindingKNM1');
%%
Time_i = M.ModelObj.TimeSec;
Time = Time_i.*(1:0.5:10);
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/results/'];
savename =[savedir,sprintf('BkgSensitivityTime_maxTime%.2f_StepSize%.2f.mat',max(Time./Time_i),(Time(2)-Time(1))./Time_i)];

if exist(savename,'file')   
    load(savename)
else
mNuSqErr = zeros(numel(Time),1);
mNuSqErrCM = zeros(numel(Time),1);
mNuSqErrNP = zeros(numel(Time),1);

for i=1:numel(Time)
    progressbar(i/numel(Time));
    M.SimulateStackRuns;
    M.ModelObj.TimeSec = Time(i);
    M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
    M.RunData.TBDIS = M.ModelObj.TBDIS;
    M.NonPoissonScaleFactor = 1.0;
    M.chi2 = 'chi2Stat';
    M.Fit;
    mNuSqErr(i) = M.FitResult.err(1);
    M.ModelObj.SetFitBias(0);
    
%     M.SimulateStackRuns;
%     M.ModelObj.TimeSec = Time(i);
%     M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
%     M.RunData.TBDIS = M.ModelObj.TBDIS;
%     M.chi2 = 'chi2Stat';
%     M.NonPoissonScaleFactor = 1.1;
%     M.Fit;
%     mNuSqErrNP(i) = M.FitResult.err(1);
%     M.ModelObj.SetFitBias(0);
    
    M.SimulateStackRuns;
    M.ModelObj.TimeSec = Time(i);
    M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
    M.RunData.TBDIS = M.ModelObj.TBDIS;
    M.chi2 = 'chi2CMShape';
    M.NonPoissonScaleFactor = 1.1;
    M.ComputeCM('SysEffects',struct('FSD','OFF'),'BkgCM','ON');
    M.Fit;
    mNuSqErrCM(i) = M.FitResult.err(1);
    M.ModelObj.SetFitBias(0);
end
save(savename,'mNuSqErrCM','mNuSqErr','Time');%,'mNuSqErrNP'
end
%% plot
Time(mNuSqErrCM<0.1) = [];
mNuSqErr(mNuSqErrCM<0.1) =[];
mNuSqErrCM(mNuSqErrCM<0.1) = [];

plotdir =  strrep(strrep(savedir,'results','plots'),'.mat','.png');

fig1 = figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
mNuSqErrSys = sqrt(mNuSqErrCM.^2-mNuSqErr.^2);
pstat = plot(Time./Time_i,mNuSqErr,'LineWidth',5);
hold on;
plot(Time./Time_i,mNuSqErrCM,'LineWidth',5);
plot(Time./Time_i+1,mNuSqErr,'--','LineWidth',3,'Color',pstat.Color);
plot(Time./Time_i,mNuSqErr+0.05,':','LineWidth',3,'Color',pstat.Color);
PrettyFigureFormat('FontSize',24);
xlabel('Time (KNM1 time)');
ylabel(sprintf('1\\sigma sensitivity on m_\\nu^2(eV^2)'));
leg = legend(sprintf('Stat'),sprintf('Stat + 10%% NP + Bkg Slope'),'Stat time shifted + 1 KNM1 time','Stat sensitivity shitfted + 0.05 eV^2');
legend boxoff
grid on
print(fig1,[plotdir,'BkgSensitivityTine.png'],'-dpng','-r450');
hold off;
xlim([min(Time./Time_i) max(Time./Time_i)]);
%%
fig2 = figure(2);
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(Time./Time_i,mNuSqErrSys./mNuSqErr,'LineWidth',3);
grid on;
PrettyFigureFormat('FontSize',24);
xlabel('Time (KNM1 time)');
ylabel(sprintf('\\sigma_{Bkg sys} / \\sigma_{Stat}'));
print(fig2,[plotdir,'BkgSensitivityTine_ratio.png'],'-dpng','-r450');
hold off;




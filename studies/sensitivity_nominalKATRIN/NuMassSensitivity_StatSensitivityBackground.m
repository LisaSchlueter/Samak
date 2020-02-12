% Calculate Sensitvity (30 months) for different background rates
% Settings
TimeSec = (124/148).*(30*24*60*60)*30;
WGTS_B_T = 0.7.*3.6;
MACE_Bmax_T =  0.7.*6;
Q_i = 18575;
range = 30;
Scan = 'OFF';
PlotFit = 'OFF'; %only used when Scan ON
chi2 = 'chi2Stat';
RecomputeFlag = 'ON';
BKG_RateSec = [0.1:0.05:1];
MACE_Ba_T = 1e-04.*7;
mNu90         = zeros(numel(BKG_RateSec),1);
TD = 'MTDcreator';
for b = 1:numel(BKG_RateSec)
    Arg = {'TimeSec',TimeSec,'Q_i',Q_i,'range',range,...
        'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,...
        'MACE_Bmax_T',MACE_Bmax_T,...
        'Scan',Scan,'chi2',chi2,...
        'TD',TD,...
        'FPD_MeanEff',0.95,'FPD_ROIlow',14,...
        'BKG_RateSec',BKG_RateSec(b)};
    S = SensitivityStudy(Arg{:});
    [~, ~, ~, ~, ~,mNu90(b),~] = S.NuMassScan;
end
save_name = sprintf('./results/StatSensitivityBackground.mat');
save(save_name,'BKG_RateSec','mNu90','MACE_Ba_T','MACE_Bmax_T','WGTS_B_T','range','Q_i','TimeSec','TD');
%% Plot
f10 = figure('Name','Sensitivity','Renderer','opengl');
set(f10, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]); %sqrt(mNu90(mNu90>0.003))*1e3

p1 = plot(BKG_RateSec(mNu90>0.003)*1e3,sqrt(mNu90(mNu90>0.003))*1e3,'LineWidth',4,'Color',rgb('CadetBlue'));
PrettyFigureFormat;
xlabel('background (mcps)')
ylabel('neutrino mass sensitivity 90% C.L. (meV)');

grid on;
set(gca,'FontSize',18);
%% save
if ~exist('./plots/png','dir')
    mkdir ./plots/png/
    mkdir ./plots/pdf/
end
print(f10,'./plots/png/StatSensitivityBKG.png','-dpng');
publish_figurePDF(f10,'./plots/pdf/StatSensitivityBKG.pdf');

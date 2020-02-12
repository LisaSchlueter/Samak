%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and Plot:
% Statstical Sensitivty as function of time
% for different background levels
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
TimeMonth = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 0.9 1 1.1 1.3 1.5 1.7 1.9 2 5 10 15 20 25 30 35 40 50 60];%2.2 2.5 3 3.2 3.5 4 4.2 4.5 5 5.5 6 7 8 9 10 12 14 16 18 20 24 26 28 30 31 32 33 35 37 40 42 44 46 48 50 52 55 58 60];
TimeSec = (30*24*60*60).*TimeMonth.*(124/148);
WGTS_B_T = [1, 0.7, 0.7, 0.7].*3.6;
MACE_Bmax_T = [1, 0.7, 0.7, 0.7].*6;
Q_i = 18575;
range = 30;
Scan = 'OFF';
PlotFit = 'OFF'; %only used when Scan ON
chi2 = 'chi2Stat';
RecomputeFlag = 'OFF';
BKG_RateSec = [0.01,0.201,0.364,1];
MACE_Ba_T = 1e-04.*[3,7, 7, 7];
MTD = {'DR30','MTDcreator','MTDcreator','MTDcreator'};
mNu90         = zeros(numel(TimeSec),numel(BKG_RateSec));

save_name = sprintf('./results/StatSensitivity_BKG_%.0f-%.0f_maxTime%.0fMonth.mat',min(BKG_RateSec)*1e3,max(BKG_RateSec)*1e3,max(TimeMonth));
if exist(save_name,'file') && strcmp(RecomputeFlag,'OFF')
    load(save_name)
else
    for b = 1:numel(BKG_RateSec)
        Arg = {'TimeSec',TimeSec(b),'Q_i',Q_i,'range',range,...
            'MACE_Ba_T',MACE_Ba_T(b),'WGTS_B_T',WGTS_B_T(b),...
            'MACE_Bmax_T',MACE_Bmax_T(b),...
            'Scan',Scan,'chi2',chi2,...
            'TD',MTD{b},...
            'FPD_MeanEff',0.95,'FPD_ROIlow',14,...
            'BKG_RateSec',BKG_RateSec(b)};
        S = SensitivityStudy(Arg{:});
        %% Loop over RunTimes
        for i=1:numel(TimeSec)
            S.TimeSec = TimeSec(i);
            S.InitializeModels; % Init Models with new RunTime
            [~, ~, ~, ~, ~,mNu90(i,b),~] = S.NuMassScan;
        end
    end
    save(save_name,'mNu90','TimeSec','BKG_RateSec','MTD','MACE_Ba_T','MACE_Bmax_T','WGTS_B_T','Q_i');
end
%% Plot
PlotTime = TimeSec.*(148/124);
leg_labels = cell(numel(BKG_RateSec,1));
%plotmNu = interp1(PlotTime(mNu90>0.003)/(30*24*60*60), sqrt(mNu90(mNu90>0.003))*1e3,TimeSec(mNu90>0.003)/(30*24*60*60),'spline');
plotmNu = NaN.*zeros(numel(TimeSec),numel(BKG_RateSec));
for b=1:numel(BKG_RateSec)
plotmNu(mNu90(:,b)>0.003,b) = sqrt(mNu90(mNu90(:,b)>0.003,b))*1e3;
leg_labels{b} = sprintf('%.0f mcps',BKG_RateSec(b)*1e3);
end

f1 = figure('Name','Sensitivity','Renderer','opengl');
set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]); %sqrt(mNu90(mNu90>0.003))*1e3

p1 = plot(PlotTime(~isnan(plotmNu(:,1)))/(30*24*60*60), plotmNu(~isnan(plotmNu(:,1)),1),'LineWidth',4,'Color',rgb('IndianRed'));
hold on;
p2 = plot(PlotTime(~isnan(plotmNu(:,2)))/(30*24*60*60), plotmNu(~isnan(plotmNu(:,2)),2),'LineWidth',4,'Color',rgb('GoldenRod'));
p3 = plot(PlotTime(~isnan(plotmNu(:,3)))/(30*24*60*60), plotmNu(~isnan(plotmNu(:,3)),3),'LineWidth',4,'Color',rgb('CadetBlue'));
p4 = plot(PlotTime(~isnan(plotmNu(:,4)))/(30*24*60*60), plotmNu(~isnan(plotmNu(:,4)),4),'LineWidth',4,'Color',rgb('RoyalBlue'));
PrettyFigureFormat;
xlabel('time (month)')
xlim([0, 60]);
xticks((0:1:60));

myXticks = strings(numel(0:1:30),1);
for i=1:numel(myXticks)
    if mod(i,2)==1 %odd
    myXticks(i) = string(i-1);
    else
        myXticks(i) = ' ';
    end
end

leg = legend([p1,p2,p3,p4],leg_labels{:});
leg.Title.String = 'background';
leg.FontSize = 18;
legend boxoff
%xticklabels(myXticks);
s = strings(5,1);
s1 = string(0:6);
set(gca,'XScale','log');
xticklabels({s1{:},s{:},'12',s{:},'18',s{:},'24',s{:},'30',s{:},s{:},s{:},s{:},s{:},'','','','','60'})
ylim([100 1650]);
%yticks(200:200:1600);
ylabel('neutrino mass sensitivity 90% C.L. (meV)');
grid on;
set(gca,'FontSize',18);
%% save
if ~exist('./plots/png','dir')
    mkdir ./plots/png/
    mkdir ./plots/pdf/
end
print(f1,'./plots/png/StatSensitivityTimeEvolution_4BKG.png','-dpng','-r350');
publish_figurePDF(f1,'./plots/png/StatSensitivityTimeEvolution_4BKG.pdf');

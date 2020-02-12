TimeSecs = (24*60*60)*(124/148).*[1 2 3 4 5 7 10 15 20 25 30 35 40 45 50 60 70 80 90 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 430 460 500 530 560 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500 1600 1700 1800 1900 1900 2000 2100 2200 2300];
%TimeSecs = (24*60*60)*(124/148).*[1 10 20 50 100 500  1000 1500 2000 2500];
MACE_Ba_T = 7*1e-04;
WGTS_B_T = 0.7*3.6;
MACE_Bmax_T = 0.7*6;
Q_i = 18575;
range = [30,60];
Scan = 'OFF';
PlotFit = 'OFF'; %only used when Scan ON
SysBudget = '09'; % 6=42 day, 07= 300 days, 09 = 900 days
RecomputeFlag = 'OFF';
ELossFlag = 'Abdurashitov';
FPD_ROIlow = 14;
FPD_MeanEff = 0.95;
BKG_RateSec = ''; % if none given: computed according to ROI
AnchorBkg6G = '';
TD = 'MTDcreator';
%% Common Plot Options

Arg = {'TimeSec',TimeSecs(1),'Q_i',Q_i,'range',range(2),...
    'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'Scan',Scan,'PlotFit',PlotFit,'SysBudget',SysBudget,...
    'RecomputeFlag',RecomputeFlag,'TD',TD,'mNuStop',2,...
    'FPD_MeanEff',FPD_MeanEff,...
    'ELossFlag',ELossFlag,'FPD_ROIlow',FPD_ROIlow,...
    'BKG_RateSec',BKG_RateSec,'AnchorBkg6G',AnchorBkg6G};

S = SensitivityStudy(Arg{:});
%%
mNu90Stat60 = zeros(numel(TimeSecs),1);
mNu90CM60   = zeros(numel(TimeSecs),1);
for i=1:numel(TimeSecs)    
    progressbar(i/numel(TimeSecs))
  S.chi2 = 'chi2Stat';
  S.TimeSec = TimeSecs(i);
  S.InitializeModels;
   [~, ~, ~, ~, ~,mNu90Stat60(i),~] = S.NuMassScan;
   
   S.chi2 = 'chi2CM';
   S.SysEffect = 'all'; S.ComputeCM;
   [~, ~, ~, ~, ~,mNu90CM60(i),~] = S.NuMassScan; 
end

%%
mNu90Stat30 = zeros(numel(TimeSecs),1);
mNu90CM30   = zeros(numel(TimeSecs),1);
S.range = 30; S.SetTD; S.InitializeModels;
for i=1:numel(TimeSecs)    
    progressbar(i/numel(TimeSecs))
  S.chi2 = 'chi2Stat';
  S.TimeSec = TimeSecs(i);
  S.InitializeModels;
   [~, ~, ~, ~, ~,mNu90Stat30(i),~] = S.NuMassScan;
   
   S.chi2 = 'chi2CM';
   S.SysEffect = 'all'; S.ComputeCM;
   [~, ~, ~, ~, ~,mNu90CM30(i),~] = S.NuMassScan; 
end
%%
f1 = figure('Renderer','opengl');
mNu90Sys60 = sqrt(mNu90CM60.^2-mNu90Stat60.^2);
pStat = plot(TimeSecs(mNu90Stat60>0.01)./((24*60*60*30)*(124/148)),sqrt(mNu90Stat60(mNu90Stat60>0.01)),'LineWidth',3);
hold on;
pSys = plot(TimeSecs(mNu90Stat60>0.01)./((24*60*60*30)*(124/148)),sqrt(mNu90Sys60(mNu90Stat60>0.01)),'LineWidth',3);
pAll =  plot(TimeSecs(mNu90Stat60>0.01)./((24*60*60*30)*(124/148)),sqrt(mNu90CM60(mNu90Stat60>0.01)),'LineWidth',3);
leg = legend([pStat, pSys, pAll],'stat','sys','total'); legend boxoff
PrettyFigureFormat;
hold off;
xlabel('time (month)')

%%
f2 = figure('Renderer','opengl');
mNu90Sys30 = sqrt(mNu90CM30.^2-mNu90Stat30.^2);
pStat = plot(TimeSecs./((24*60*60)*(124/148)),sqrt(mNu90Stat30),'LineWidth',3);
hold on;
pSys = plot(TimeSecs./((24*60*60)*(124/148)),sqrt(mNu90Sys30),'LineWidth',3);
pAll =  plot(TimeSecs./((24*60*60)*(124/148)),sqrt(mNu90CM30),'LineWidth',3);
leg = legend([pStat, pSys, pAll],'stat','sys','total'); legend boxoff
PrettyFigureFormat;
hold off;




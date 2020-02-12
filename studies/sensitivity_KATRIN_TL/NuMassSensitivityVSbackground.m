%% 1. Step: Initialize SensitivityStudy Object 

%% KATRIN SETTINGS
TimeDays      = 42;
TimeSec       = (TimeDays*24*60*60)*(124/148);
MACE_Ba_T     = 7*1e-04;        % B-field AP
WGTS_B_T      = 0.7*3.6;        % B-field WGTS
MACE_Bmax_T   = 0.7*6; 
Q_i           = 18575;          % Endpoint
range         = 30;             % Max range, V
TD            = 'MTDcreator';   % Type of MTD
Anchor6GValue = 500;         % Background @6G in AP

%% COMPUTATION SETTINGS
RecomputeFlag = 'ON';
Scan          = 'ON';
ScanPrcsn     = 0.02;
mNuStop       = 2;                
PlotFit       = 'OFF';  

% SYSTEMATICS
SysBudget     = '06';
chi2          = 'chi2CMShape';
SysEffect     = 'all';
InputArg = {...
    'TimeSec',TimeSec,...
    'Q_i',Q_i,...
    'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'SysEffect',SysEffect,...
    'SysBudget',SysBudget,...
    'TD',TD,...
    'range',range,...
    'Scan',Scan,...
    'ScanPrcsn',ScanPrcsn,...
    'mNuStop',mNuStop,...
    'RecomputeFlag',RecomputeFlag...
    'PlotFit',PlotFit,...
    };
%S=SensitivityStudy(InputArg{:},'AnchorBkg6G',450*1e-3);S.NuMassScan;return;

% CREATE AND LOOP ON SENSITIVITY OBJECT
bkgV = 10:100:1000;
mNu90_stat = zeros(1,numel(bkgV));
mNu90_sys = zeros(1,numel(bkgV));
counter=0;
for b=bkgV
    counter=counter+1;
    Sstat=SensitivityStudy(InputArg{:},'AnchorBkg6G',b,'chi2','chi2Stat');
    Ssys=SensitivityStudy(InputArg{:},'AnchorBkg6G',b,'chi2','chi2CM');
    fprintf(2,'******************  Computing with Background = %.0f mcps ******************\n',b);
    [~, ~, ~, ~,~,mNu90_stat(counter),~] = Sstat.NuMassScan;    
    [~, ~, ~, ~,~,mNu90_sys(counter),~] = Ssys.NuMassScan;    
    fprintf(2,'****************************************************************************\n');
    fprintf(2,'*************  Background = %.0f mcps  - mNu90 = %.3f     ******************\n',b,mNu90_stat(counter));
    fprintf(2,'****************************************************************************\n');
end
disp(mNu90_stat)
disp(mNu90_sys)

%%
plt = Plot(bkgV,smooth(sqrt(mNu90_stat))*1e3,bkgV,smooth(sqrt(mNu90_sys))*1e3);
plt.Colors = {rgb('DarkBlue'),rgb('CadetBlue')};       
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'stat only',sprintf('stat+sys (%s)',SysBudget)}; % legends
plt.LegendLoc  =  'NorthWest';
plt.XLabel = 'Background (mcps)'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
mytitle = sprintf('KATRIN - %.0f days - MTD -%0.f V - Eff_{FPD} x Cov_{FPD} = %0.2f , B_A=%0.1f G , B_s=%0.1f T , B_m=%0.1f T (Samak)',TimeDays,range,FPD_MeanEff,MACE_Ba_T*1e4,WGTS_B_T,MACE_Bmax_T);
plt.Title  = mytitle; % plot title
plt.LineWidth = [2, 5]; % three line widths
plt.LineStyle = {'--', '-'}; % three line styles
plt.MarkerSpacing = [10, 10];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 20;
outputfile = sprintf('KATRIN_Background_%.0fdays_MTD%0.f_EffCovFPD%0.2fBA%0.1fGBs%0.1fTBm%0.1fT',TimeDays,range,FPD_MeanEff,MACE_Ba_T*1e4,WGTS_B_T,MACE_Bmax_T);
plt.export([outputfile '.png']);
save([outputfile '.mat'],'mNu90_stat','mNu90_sys','TimeDays','TimeSec','range','FPD_MeanEff','MACE_Ba_T','WGTS_B_T','MACE_Bmax_T','SysBudget','Q_i','FPD_MeanEff','Sstat','Ssys');

%
%
%
%

% Common Option
commonOpt = {'TD','Sensitivity_60eV_Ba6G','WGTS_B_T',3.6*0.7,'MACE_Ba_T',6e-04,'mnuSq_i_Fit',(0:21:1000)*1e-03,'plotFit','OFF','SysEffect','all',};

% %% TDR 2004 - 5 years plan
% % Stat 
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2019_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',219*86400,commonOpt{:},'chi2','chi2Stat');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2020_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',2*219*86400,commonOpt{:},'chi2','chi2Stat');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2021_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',3*219*86400,commonOpt{:},'chi2','chi2Stat');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2022_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',4*219*86400,commonOpt{:},'chi2','chi2Stat');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2023_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',5*219*86400,commonOpt{:},'chi2','chi2Stat');
% mNu90_tdr_s = [mNu90_tdr_2019_s mNu90_tdr_2020_s mNu90_tdr_2021_s mNu90_tdr_2022_s mNu90_tdr_2023_s];
% 
% % Stat+Syst
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2019_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',219*86400,commonOpt{:},'chi2','chi2CM');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2020_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',2*219*86400,commonOpt{:},'chi2','chi2CM');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2021_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',3*219*86400,commonOpt{:},'chi2','chi2CM');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2022_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',4*219*86400,commonOpt{:},'chi2','chi2CM');
% [mnuSq_i_Fit, chi2minScan, dofScan,mNu90_tdr_2023_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',5*219*86400,commonOpt{:},'chi2','chi2CM');
% mNu90_tdr_ss = [mNu90_tdr_2019_ss mNu90_tdr_2020_ss mNu90_tdr_2021_ss mNu90_tdr_2022_ss mNu90_tdr_2023_ss];

%% Nominal - 5 years plan
% Stat 
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2019_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',179.5*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2020_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',2*179.5*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2021_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',3*179.5*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2022_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',4*179.5*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2023_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',5*179.5*86400,commonOpt{:},'chi2','chi2Stat');
mNu90_nom_s = [mNu90_nom_2019_s mNu90_nom_2020_s mNu90_nom_2021_s mNu90_nom_2022_s mNu90_nom_2023_s];

% Stat+Syst
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2019_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',179.5*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2020_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',2*179.5*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2021_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',3*179.5*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2022_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',4*179.5*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_nom_2023_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',5*179.5*86400,commonOpt{:},'chi2','chi2CM');
mNu90_nom_ss = [mNu90_nom_2019_ss mNu90_nom_2020_ss mNu90_nom_2021_ss mNu90_nom_2022_ss mNu90_nom_2023_ss];

%% Reduced - 5 years plan
% Stat 
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2019_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',42*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2020_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+1*102)*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2021_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+2*102)*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2022_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+3*102)*86400,commonOpt{:},'chi2','chi2Stat');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2023_s,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+4*102)*86400,commonOpt{:},'chi2','chi2Stat');
mNu90_red_s = [mNu90_red_2019_s mNu90_red_2020_s mNu90_red_2021_s mNu90_red_2022_s mNu90_red_2023_s];

% Stat+Syst
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2019_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',42*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2020_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+1*102)*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2021_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+2*102)*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2022_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+3*102)*86400,commonOpt{:},'chi2','chi2CM');
[mnuSq_i_Fit, chi2minScan, dofScan,mNu90_red_2023_ss,mNumin] = SensitivityNominal_NuMassScan('TimeSec',(42+4*102)*86400,commonOpt{:},'chi2','chi2CM');
mNu90_red_ss = [mNu90_red_2019_ss mNu90_red_2020_ss mNu90_red_2021_ss mNu90_red_2022_ss mNu90_red_2023_ss];

years         = [2019 2020 2021 2022 2023];
mNu90_nom_s   = sqrt(mNu90_nom_s);
mNu90_nom_ss  = sqrt(mNu90_nom_ss);

mNu90_red_s   = sqrt(mNu90_red_s);
mNu90_red_ss  = sqrt(mNu90_red_ss);
 
%% Plot Results 60_eV

%  mNu90_tdr_s   = sqrt([ 0.11055    0.063631    0.044427    0.034094    0.027652 ]);
%  mNu90_tdr_ss  = sqrt([ 0.26929     0.20061     0.18593     0.17465     0.16476 ]);
%  
%  mNu90_nom_s   = sqrt([ 0.12954    0.077465    0.054766    0.042288    0.034423 ]);
%  mNu90_nom_ss  = sqrt([ 0.26929     0.21189     0.19417     0.18393     0.17509 ]);
%  
%  mNu90_red_s   = sqrt([ 0.2376     0.14492     0.10153    0.077245    0.062164 ]);
%  mNu90_red_ss  = sqrt([ 0.43683     0.29047     0.23556      0.2117     0.19952 ]);
 
%  mNu90_tdr_s   = sqrt(mNu90_tdr_s);
%  mNu90_tdr_ss  = sqrt(mNu90_tdr_ss);
 
% Plot Results 

 

% Nominal
figure(1)
plot(years,mNu90_nom_s,years,mNu90_nom_ss);
plt = Plot();
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'statistical limit', 'stat+sys (samak)'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN (Nominal 5 x 180 days) - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [2, 4]; % three line widths
plt.LineStyle = {'--', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.1 1];

% TDR
% figure(2)
% plot(years,mNu90_tdr_s,years,mNu90_tdr_ss);
% plt = Plot();
% plt.YGrid = 'on'; % 'on' or 'off'
% plt.BoxDim = [14, 6]; %[width, height] in inches
% plt.Legend = {'statistical limit', 'stat+sys (samak)'}; % legends
% plt.XLabel = 'Year'; % xlabel
% plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
% plt.Title  = 'KATRIN (TDR 5 x 219 days) - B_A=6 G \rightarrow B=230 mcps - 60 eV Optimized MTD'; % plot title
% plt.LineWidth = [2, 4]; % three line widths
% plt.LineStyle = {'--', '-'}; % three line styles
% plt.MarkerSpacing = [10, 100];
% plt.YScale = 'lin'; % 'linear' or 'log'
% plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
% plt.YLim  = [0.1 0.7];

% Reduced
figure(3)
plot(years,mNu90_red_s,years,mNu90_red_ss);
plt = Plot();
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'statistical limit', 'stat+sys (samak)'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN (Reduced 42 + 4 x 102 days) - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [2, 4]; % three line widths
plt.LineStyle = {'--', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.1 1];

% stat
figure(4)
plot(years,mNu90_nom_s,years,mNu90_red_s);
plt = Plot();
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'Nominal Stat', 'Reduced Stat'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN statistical limit - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [4, 4]; % three line widths
plt.LineStyle = {'-', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.1 1];

% sys Samak
figure(5)
plot(years,mNu90_nom_s,years,mNu90_red_s,years,mNu90_nom_ss,years,mNu90_red_ss);
plt = Plot();
plt.Colors = { % three colors for three data set
    [1, 0, 0] % data set 1
    [0.25, 0.25, 0.25] % data set 2
    [1, 0, 0] % data set 3
    [0.25, 0.25, 0.25]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'Nominal Stat - stat only', 'Reduced Stat - stat only','Nominal Stat - stat + sys', 'Reduced Stat - stat + sys'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'm_\nu Upper Limit - 90% CL'; %ylabel
plt.Title  = 'KATRIN stat + sys (Samak) - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [2, 2, 4, 4]; % three line widths
plt.LineStyle = {'--', '--','-', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.1 1];
plt.export('Katrin60eVMTDNominalReducedStat_SamakSys.png');

% sys Samak
figure(55)
sys1sm2= 0.017;
stat1sm2_nom = ((mNu90_nom_s)./(1.64).^0.5).^2;
stat1sm2_red = ((mNu90_red_s)./(1.64).^0.5).^2;
mNu90_nom_sstdr = (1.64.*(stat1sm2_nom.^2+sys1sm2.^2).^0.5).^0.5;
mNu90_red_sstdr = (1.64.*(stat1sm2_red.^2+sys1sm2.^2).^0.5).^0.5;
plot(years,mNu90_nom_s,years,mNu90_red_s,years,mNu90_nom_sstdr,years,mNu90_red_sstdr);
plt = Plot();
plt.Colors = { % three colors for three data set
    [1, 0, 0] % data set 1
    [0.25, 0.25, 0.25] % data set 2
    [1, 0, 0] % data set 3
    [0.25, 0.25, 0.25]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {...
    'Nominal statistics (5 x 180 days) - stat only', ...
    'Reduced statistics (42 + 4 x 102 days) - stat only',...
    'Nominal - stat + sys (TDR2005) ', ...
    'Reduced - stat + sys (TDR2005)'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'm_\nu Sensitivity - 90% CL'; %ylabel
plt.Title  = 'KATRIN stat + sys (TDR2005) - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [2, 2, 4, 4]; % three line widths
plt.LineStyle = {'--', '--','-', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.2 0.6];
plt.export('Katrin60eVMTDNominalReducedStatTDRSys.png');


% sys / relative to final sensitivity
figure(6)
plot(years,mNu90_nom_s./mNu90_nom_s(end)*1,...
    years,mNu90_red_s./mNu90_red_s(end)*1,...
    years,mNu90_nom_ss./mNu90_nom_ss(end)*1,...
    years,mNu90_red_ss./mNu90_red_ss(end)*1);
plt = Plot();
plt.Colors = { % three colors for three data set
    [1, 0, 0] % data set 1
    [0.25, 0.25, 0.25] % data set 2
    [1, 0, 0] % data set 3
    [0.25, 0.25, 0.25]};
plt.YGrid = 'on'; % 'on' or 'off'
plt.BoxDim = [14, 6]; %[width, height] in inches
plt.Legend = {'m_\nu Upper Limit - 90% CL', 'Reduced Stat - stat only','Nominal Stat - stat + sys', 'Reduced Stat - stat + sys'}; % legends
plt.XLabel = 'Year'; % xlabel
plt.YLabel = 'Sensitivity / Final Sensitivity'; %ylabel
plt.Title  = 'KATRIN stat + sys (Samak) - B_S = 2.5 T - B_A=6 G \rightarrow B=230 mcps - 60 eV Opt. MTD'; % plot title
plt.LineWidth = [2, 2, 4, 4]; % three line widths
plt.LineStyle = {'--', '--','-', '-'}; % three line styles
plt.MarkerSpacing = [10, 100];
plt.YScale = 'lin'; % 'linear' or 'log'
plt.XTick = [2019,2020,2021,2022,2023]; %[tick1, tick2, .. ]
plt.YLim  = [0.9 2];
plt.export('Katrin60eVMTDNominalReducedStat-Relative.png');

WGTS_CD_MolPerCm2 = 4.4531e+17;
run = 40680;
%run = 11111;

% Build Model
S100 = ref_runsummaries(run,'ex2','ISCS','Theory','recomputeRF','OFF','WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2);
S100.TimeSec = 3*3600;
S100.ComputeTBDDS;S100.ComputeTBDIS

% Build Model
S75 = ref_runsummaries(run,'ex2','ISCS','Theory','recomputeRF','OFF','WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2*0.75);
S75.TimeSec = 4.5*3600;
S75.ComputeTBDDS;S75.ComputeTBDIS

% Build Model
S50 = ref_runsummaries(run,'ex2','ISCS','Theory','recomputeRF','OFF','WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2*0.5);
S50.TimeSec = 6*3600;
S50.ComputeTBDDS;S50.ComputeTBDIS

% Build Model
S25 = ref_runsummaries(run,'ex2','ISCS','Theory','recomputeRF','OFF','WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2*0.25);
S25.TimeSec = 12*3600;
S25.ComputeTBDDS;S25.ComputeTBDIS

%% Plot Counts
clear plt;
plt = Plot(S100.qU-S100.Q,S100.TBDIS, S100.qU-S100.Q,S75.TBDIS, S100.qU-S100.Q,S50.TBDIS, S100.qU-S100.Q,S25.TBDIS);
plt.Title = 'Integral Spectra : Colmun Density Study - 22 hours'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Counts'; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
plt.XLim   = [-1650 50];
%plt.XLim   = [-10000 50];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 4.4e17 cm^{-2}', '\rho d = 3.3e17 cm^{-2}', '\rho d = 2.2e17 cm^{-2}','\rho d = 1.1e17 cm^{-2}'}; % legends
plt.LegendLoc = 'SouthWest';
plt.export('rhDdatataking-wholeRange-Count.png');
PrettyFigureFormat

% Plot Rates
clear plt;
plt = Plot(S100.qU-S100.Q,S100.TBDIS./S100.qUfrac./S100.TimeSec,...
 S100.qU-S100.Q,S75.TBDIS./S75.qUfrac./S75.TimeSec,...
 S100.qU-S100.Q,S50.TBDIS./S50.qUfrac./S50.TimeSec,...
 S100.qU-S100.Q,S25.TBDIS./S25.qUfrac./S25.TimeSec);
plt.Title = 'Integral Spectra : Colmun Density Study - 22 hours'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Rate'; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
plt.XLim   = [-1650 50];
%plt.XLim   = [-10000 50];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 4.4e17 cm^{-2}', '\rho d = 3.3e17 cm^{-2}', '\rho d = 2.2e17 cm^{-2}','\rho d = 1.1e17 cm^{-2}'}; % legends
plt.LegendLoc = 'SouthWest';
plt.export('rhDdatataking-wholeRange-Count.png');
PrettyFigureFormat

%% Plot Ratio to S100
clear plt;
plt = Plot(S100.qU-S100.Q,(S100.TBDIS./S100.qUfrac./S100.TimeSec)./(S100.TBDIS./S100.qUfrac./S100.TimeSec),...
 S100.qU-S100.Q,(S75.TBDIS./S75.qUfrac./S75.TimeSec)./(S100.TBDIS./S100.qUfrac./S100.TimeSec),...
 S100.qU-S100.Q,(S50.TBDIS./S50.qUfrac./S50.TimeSec)./(S100.TBDIS./S100.qUfrac./S100.TimeSec),...
 S100.qU-S100.Q,(S25.TBDIS./S25.qUfrac./S25.TimeSec)./(S100.TBDIS./S100.qUfrac./S100.TimeSec));
plt.Title = 'Integral Spectra : Colmun Density Study - 22 hours'; % plot title
plt.XLabel = 'E-qU (V)'; % xlabel
plt.YLabel = 'Rate'; %ylabel
plt.YScale = 'log'; % 'linear' or 'log'
plt.XLim   = [-1650 50];
%plt.XLim   = [-10000 50];
plt.XScale = 'lin'; % 'linear' or 'log'
plt.FontSize = 16;
plt.Legend = {'\rho d = 4.4e17 cm^{-2}', '\rho d = 3.3e17 cm^{-2}', '\rho d = 2.2e17 cm^{-2}','\rho d = 1.1e17 cm^{-2}'}; % legends
plt.LegendLoc = 'SouthWest';
plt.export('rhDdatataking-wholeRange-Ratio.png');
PrettyFigureFormat


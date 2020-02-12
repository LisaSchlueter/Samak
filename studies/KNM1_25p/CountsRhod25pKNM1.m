% Generic Options
% Data Taking Time
TimeDays = 30; TimeSec = TimeDays*24*60*60;

% Pixel List for Un iform Fit
FPD_SegOffPixList = 1:119;
FPD_ROIlow        = 14;
FPD_MeanEff       = 0.95;

% Column Density
WGTS_CD_MolPerCm2_22  = 1.1e17;
WGTS_CD_MolPerCm2_30  = 5e17*0.3;
WGTS_CD_MolPerCm2_50  = 2.5e17;
WGTS_CD_MolPerCm2_100 = 5e17;

% ELoss
ELossFlag         = 'KatrinD2';

% Magnetic Fields
MACE_Ba_T         = 6*1e-04;
WGTS_B_T          = 0.7*3.6;
MACE_Bmax_T       = 0.7*6;

% EndPoint
Q_i               = 18574;

% MTD
TD = 'KNM1_25p';

% Background
                BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,...
                    'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',FPD_ROIlow,...
                    'AnchorBkg6G',430e-3, ...
                    'FPD_SegOffPixList',FPD_SegOffPixList);
 
ModelArg = {...
    'FPD_SegOffPixList',FPD_SegOffPixList,...
    'BKG_RateAllFPDSec',BKG_RateSec,...
    'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeSec',TimeSec,...
    'TD',TD,...
    'mnuSq_i',0,...
    'Q_i',Q_i,...
    'FPD_MeanEff',FPD_MeanEff,...
    'ELossFlag',ELossFlag,...
    'recomputeRF','ON',...
    'FPD_ROIlow',FPD_ROIlow};
            
SimObj22  = ref_TBD_NominalKATRIN(ModelArg{:},'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_22);
SimObj30  = ref_TBD_NominalKATRIN(ModelArg{:},'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_30);
SimObj50  = ref_TBD_NominalKATRIN(ModelArg{:},'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_50);
SimObj100 = ref_TBD_NominalKATRIN(ModelArg{:},'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_100);
%%
SimObj22.ComputeTBDDS;  SimObj22.ComputeTBDIS;
SimObj30.ComputeTBDDS;  SimObj30.ComputeTBDIS;
SimObj50.ComputeTBDDS;  SimObj50.ComputeTBDIS;
SimObj100.ComputeTBDDS; SimObj100.ComputeTBDIS;


%%
p1 = Plot(SimObj22.qU(2:end)-Q_i,...
    cumsum(SimObj22.TBDIS(2:end),'reverse'),...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj50.TBDIS(2:end),'reverse'),...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj100.TBDIS(2:end),'reverse'));
p1.YScale = 'log';
p1.BoxDim = [10 , 5];

%%
p2 = Plot(SimObj22.qU(2:end)-Q_i,...
    cumsum(SimObj22.TBDIS(2:end),'reverse')./cumsum(SimObj22.TBDIS(2:end),'reverse')*100,...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj50.TBDIS(2:end),'reverse')./cumsum(SimObj100.TBDIS(2:end),'reverse')*100,...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj22.TBDIS(2:end),'reverse')./cumsum(SimObj100.TBDIS(2:end),'reverse')*100);
p2.YScale = 'log';
p2.BoxDim = [10 , 5];
p2.Colors = {[0, 0, 0],rgb('SteelBlue'),rgb('Amethyst')}; % [red, green, blue]
p2.LineWidth = 3; % line width
p2.LineStyle = {'--','-','-'}; % line style: '-', ':', '--' etc
p2.Legend = {'100% \rho d', '50% \rho d', '22% \rho d'}; % legends
            ylabel('rate (cps)','FontSize',16);
            p2.XLabel = 'retarding potential - 18574 (V)';  
p2.YLabel = 'Cumulative Counts in Percent';  
p2.XLim = [-30, 10]; % [min, max]
p2.Title = 'Cumulative Statistics for KNM1';
p2.export('CountsRhod25pKNM1_1.png');

%%
p3 = Plot(SimObj22.qU(2:end)-Q_i,...
    cumsum(SimObj22.TBDIS(2:end),'reverse')./cumsum(SimObj22.TBDIS(2:end),'reverse')*100,...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj30.TBDIS(2:end),'reverse')./cumsum(SimObj100.TBDIS(2:end),'reverse')*100,...
    SimObj100.qU(2:end)-Q_i,...
    cumsum(SimObj22.TBDIS(2:end),'reverse')./cumsum(SimObj100.TBDIS(2:end),'reverse')*100);
p3.YScale = 'log';
p3.BoxDim = [10 , 5];
p3.Colors = {[0, 0, 0],rgb('SteelBlue'),rgb('Amethyst')}; % [red, green, blue]
p3.LineWidth = 3; % line width
p3.LineStyle = {'--','-','-'}; % line style: '-', ':', '--' etc
p3.Legend = {'100% \rho d', '30% \rho d', '22% \rho d'}; % legends
            ylabel('rate (cps)','FontSize',16);
            p3.XLabel = 'retarding potential - 18574 (V)';  
p3.YLabel = 'Cumulative Counts in Percent';  
p3.XLim = [-30, 10]; % [min, max]
p3.Title = 'Cumulative Statistics for KNM1';
p3.export('CountsRhod25pKNM1_2.png');

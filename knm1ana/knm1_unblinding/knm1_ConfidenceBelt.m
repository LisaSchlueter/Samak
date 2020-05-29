% KNM1 confidence belt construction
% calculated or imported from data
% calculation and plotting takes place in RunSensitivity class
% Lisa Schlüter, 2019
%% settings
Mode        = 'FC';  % FC = Feldman Cousin, LT = Lokov Tkachov
Sensitivity = 'OFF'; % OFF= show best fit, ON = show sensitivity only
SavePlot    = 'ON';

RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...%free par
    'exclDataStart',13,...
    'SysBudget',22,...
    'FSDFlag','Sibille0p5eV',...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1.064,...
    'RadiativeFlag','ON',...
    'ELossFlag','KatrinT2',...
    'chi2','chi2CMShape'};

% set up model 
D = MultiRunAnalysis(RunAnaArg{:});
S = RunSensitivity('RunAnaObj',D);
%% Define Test values
mNuSq_1 = [0,0.15,0.2:0.1:0.7,0.9,1.2,1.35,1.5];
mNuSq_2 = [0.11,0.13,0.17,0.25:0.1:0.95,1,1.1,1.3,1.6,1.7,1.8];
mNuSq_3 = [0.8,1.4];
mNuSq = sort([mNuSq_1,mNuSq_2,mNuSq_3]);
%% Compute and Plot Confidence Belt
S.ConfLevel = 0.9; % confidence level (0==1 sigma)
switch Mode
    case 'FC'
        S.ComputeFC_Asimov('mNuSq_t',mNuSq,'nSamplesAsimov',300);
        S.PlotFCBelt('HoldOn','OFF','Sensitivity',Sensitivity,'SavePlot',SavePlot);
    case 'LT'
        S.ComputeLokhov_Asimov('mNuSq_t',mNuSq);
        S.PlotFCBelt('Lokov','ON','Sensitivity',Sensitivity,'SavePlot',SavePlot);
end
%% plot likelihood function
S.PlotFC_DeltaChi2('PDF','Central','SavePlot','ON','mNuSq_t',0); % probability density function with best fit probability
S.PlotFC_DeltaChi2('PDF','1sigma','SavePlot','ON','mNuSq_t',0);  % probability density function with 1 sigma boundaries
S.PlotFC_DeltaChi2('PDF','OFF','SavePlot','ON','mNuSq_t',0.2);   % Delta chi2 curve

%S.PlotFC_PDF('KNM1Central','OFF','mNuSq_t',0);


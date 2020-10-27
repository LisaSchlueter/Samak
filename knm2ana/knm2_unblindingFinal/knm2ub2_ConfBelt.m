% KNM1 confidence belt construction
% calculated or imported from data
% calculation and plotting takes place in RunSensitivity class
% Lisa Schlüter, 2019
%% settings
Mode        = 'FC';  % FC = Feldman Cousin, LT = Lokov Tkachov
Sensitivity = 'ON'; % OFF= show best fit, ON = show sensitivity only
SavePlot    = 'ON';
range = 40;
chi2 = 'chi2CMShape';
SysBudget = 38;
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
SigmaSq =  0.0124+0.0025;
RunAnaArg = {'RunList','KNM2_Prompt',...
    'fixPar','mNu E0 Norm Bkg',...%free par
    'SysBudget',SysBudget,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(SigmaSq),...
    'TwinBias_FSDSigma',sqrt(SigmaSq),...
    'TwinBias_Q',18573.7};

%% set up model
D = MultiRunAnalysis(RunAnaArg{:});
D.exclDataStart = D.GetexclDataStart(range);
S = RunSensitivity('RunAnaObj',D);
%% Define Test values
% mNuSq_1 = [0,0.05,0.08,0.15,0.2:0.1:0.9,1.2,1.35,1.5];
% mNuSq_2 = [0.11,0.13,0.17,0.25:0.1:0.55,0.95];
% mNuSq_3 = [1,1.1,1.3,1.6,1.4,1.7];
%mNuSq = sort([mNuSq_1,mNuSq_2,mNuSq_3]);
mNuSq = sort([0:0.1:1.6,0.05,0.08]);
%mNuSq = mNuSq(1:15);
%% Compute and Plot Confidence Belt
S.ConfLevel = 0.9; % confidence level (0==1 sigma)
switch Mode
    case 'FC'
        S.ComputeFC_Asimov('mNuSq_t',mNuSq,'nSamplesAsimov',300);
        S.PlotFCBelt('HoldOn','OFF','Sensitivity','OFF',...%Sensitivity,...
            'SavePlot','ON','XLim',[-1.5,1.5],...
            'Style','Pretty','mNuSq_bf',0.190);
    case 'LT'
        mNuSq(mNuSq==0.8) = [];
        S.ComputeLokhov_Asimov('mNuSq_t',mNuSq);
        S.PlotFCBelt('Lokov','ON','Sensitivity','OFF','SavePlot',SavePlot,...
            'Style','Pretty','XLim',[-1.5,1.5],'mNuSq_bf',0.190);
end
%% plot likelihood function
%S.PlotFC_DeltaChi2('PDF','1sigma','SavePlot','ON','mNuSq_t',0);  % probability density function with 1 sigma boundaries

%S.PlotFC_DeltaChi2('PDF','Central','SavePlot','ON','mNuSq_t',0); % probability density function with best fit probability
%S.PlotFC_DeltaChi2('PDF','OFF','SavePlot','ON','mNuSq_t',0.2);   % Delta chi2 curve
%S.PlotFC_PDF('KNM1Central','OFF','mNuSq_t',0);


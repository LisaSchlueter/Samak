% KNM1 confidence belt construction
% calculated or imported from data
% calculation and plotting takes place in RunSensitivity class
% Re-Analysis with updated model parameters + uncerainties
% Lisa Schlüter, 2021
%% settings
Mode        = 'LT';  % FC = Feldman Cousin, LT = Lokov Tkachov
Sensitivity = 'OFF'; % OFF= show best fit, ON = show sensitivity only
SavePlot    = 'OFF';

RunAnaArg = {'RunList','KNM1',...
    'DataType','Twin',...
    'fixPar','mNu E0 Norm Bkg',...%free par
    'exclDataStart',13,...
    'SysBudget',200,...
    'FSDFlag','KNM2_0p5eV',...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1.064,...
    'RadiativeFlag','ON',...
    'SynchrotronFlag','ON',...
    'AngularTF','ON',...
    'ELossFlag','KatrinT2A20',...
    'chi2','chi2CMShape',...
    'BKG_PtSlope', -2.2*1e-06,...
    'TwinBias_BKG_PtSlope',-2.2*1e-06,...
    'DopplerEffectFlag','FSD'};

% set up model 
D = MultiRunAnalysis(RunAnaArg{:});

%% re-init again
D = MultiRunAnalysis(RunAnaArg{:});
S = RunSensitivity('RunAnaObj',D);
%% Define Test values
mNuSq = [1.3,1.32,1.33,1.35];%[1.3,1.32,1.35,1.4];
%% Compute and Plot Confidence Belt
S.ConfLevel = 0.9; % confidence level (0==1 sigma)
switch Mode
    case 'FC'
        S.ComputeFC_Asimov('mNuSq_t',mNuSq,'nSamplesAsimov',300);
        S.PlotFCBelt('HoldOn','OFF','Sensitivity',Sensitivity,'SavePlot',SavePlot);
    case 'LT'
        S.ComputeLokhov_Asimov('mNuSq_t',mNuSq);
     %  S.PlotFCBelt('Lokov','ON','Sensitivity',Sensitivity,'SavePlot',SavePlot);
end

%%
%S.FC_x1
GetFigure

% comarison to regular KNM1_
savedir = [getenv('SamakPath'),'knm1ana/knm1_unblinding/results/'];
savename = sprintf('%sknm1_ConfidenceBelt_LT.mat',savedir);
d1 = importdata(savename);
%d1.PlotFCBelt('Lokov','ON','Sensitivity',Sensitivity,'SavePlot','OFF');

pre = plot(S.FC_x1,mNuSq,'.-.','LineWidth',2,'MarkerSize',15,'Color',rgb('Orange'));
hold on 
p1 = plot(d1.FC_x1,d1.FC_mNuSqTrue,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
grid on
xlabel(sprintf('{\\itm}_\\nu^2 measured (eV^2)'));
ylabel(sprintf('{\\itm}_\\nu^2 true (eV^2)'));
xlim([-0.1 0.1])
%ylim(sqrt([1.12 1.5]))
PrettyFigureFormat;

SensLim1 = sqrt(interp1(d1.FC_x1,d1.FC_mNuSqTrue,0,'spline'));
SensLim2 = sqrt(interp1(S.FC_x1,S.FC_mNuSqTrue,0,'spline'));
leg = legend([p1,pre],sprintf('KNM-1 PRL Sensitivity Limit:      {\\itm}_\\nu < %.3f eV',SensLim1),...
                sprintf('KNM-1 Re-Ana Sensitivity Limit: {\\itm}_\\nu < %.3f eV',SensLim2),...
               'Location','southeast');
PrettyLegendFormat(leg);
%%
mNuSqIdx = 2;
ExclLog = S.FC_mNuSqFit(mNuSqIdx,:)<0;
PDF = S.FC_PDF(mNuSqIdx,:)./simpsons(S.FC_mNuSqFit(mNuSqIdx,:),S.FC_PDF(mNuSqIdx,:));

[CumProb,ia] = unique(S.FC_CumProb(mNuSqIdx,:));
interp1(S.FC_mNuSqFit(mNuSqIdx,ia),1-CumProb,0);

GetFigure
area(S.FC_mNuSqFit(mNuSqIdx,~ExclLog),PDF(~ExclLog),...
    'FaceColor',rgb('Silver'),'EdgeColor','k')
hold on;
area(S.FC_mNuSqFit(mNuSqIdx,ExclLog),PDF(ExclLog),...
    'FaceColor',rgb('White'),'EdgeColor','k');
grid on
xlabel(sprintf('{\\itm}_\\nu^2 measured (eV^2)'))
PrettyFigureFormat;
xlim([-2 5]);
leg = legend(sprintf('MC truth: {\\itm}_\\nu^2 = %.2f eV^2',mNuSq(mNuSqIdx)));
PrettyLegendFormat(leg);






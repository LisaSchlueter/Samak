% attempt to extract uncertainties on best-fit 


%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';%'Real';
nGridSteps = 50;
range = 40;
Mode = 'Compute';
LocalFontSize = 24;
freePar ='E0 Norm Bkg';
% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5}};
%%
S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
%

S.LoadGridFile(S.LoadGridArg{:},'ExtmNu4Sq','ON'); 
S.Interp1Grid('nInter',1e3);
S.FindBestFit;
mNu4Sq_bf = S.mNu4Sq_bf;
sin2T4_bf = S.sin2T4_bf;

%%
% get error on m4^2
mNu4Sq = S.mNu4Sq(1,:);
sin2T4 =  S.sin2T4(:,1);
chi2min = zeros(numel(mNu4Sq),1);

for i=1:numel(mNu4Sq)
  chi2min(i) = min(S.chi2(:,i));
end

%% get 1 sigma error
Idx1_inter = mNu4Sq>=10 & mNu4Sq<=S.mNu4Sq_bf;
Idx2_inter = mNu4Sq>S.mNu4Sq_bf & mNu4Sq<200 ;

mNu4Sq_errlow = S.mNu4Sq_bf-interp1(chi2min(Idx1_inter)-S.chi2_bf,mNu4Sq(Idx1_inter),1,'spline');
mNu4Sq_errup = interp1(chi2min(Idx2_inter)-S.chi2_bf,mNu4Sq(Idx2_inter),1,'spline')-S.mNu4Sq_bf;
%%
close all;
GetFigure;
plot(mNu4Sq,ones(numel(mNu4Sq),1),'k:','LineWidth',2);
hold on;
p1 = plot(mNu4Sq,chi2min-S.chi2_bf,'LineWidth',2);
hold off
set(gca,'XScale','log');
xlabel(sprintf('{\\itm}_4^2'));
ylabel(sprintf('\\chi^2 - \\chi^2_{min}'));
PrettyFigureFormat('FontSize',22);
xlim([0.1 1600]);
%%
leg = legend(p1,sprintf('Best fit: {\\itm}_4^2 = %.1f (-%.1f +%.1f) eV^2',S.mNu4Sq_bf,mNu4Sq_errlow,mNu4Sq_errup),...
    'Location','northwest');
PrettyLegendFormat(leg);

%print(gcf,'./plots/ksn2_m4Sq_bf_err.png','-dpng','-r450');

%% get error on Ue4^2
mNu4Sq = S.mNu4Sq(1,:);
sin2T4 =  S.sin2T4(:,1);
chi2min = zeros(numel(mNu4Sq),1);

for i=1:numel(mNu4Sq)
  chi2min(i) = min(S.chi2(i,:));
end
%%
% get 1 sigma error
Idx1_inter = sin2T4<=S.sin2T4_bf;
Idx2_inter = sin2T4>S.sin2T4_bf & sin2T4<0.1;

sin2T4_errlow = S.sin2T4_bf-interp1(chi2min(Idx1_inter)-S.chi2_bf,sin2T4(Idx1_inter),1,'spline');
sin2T4_errup = interp1(chi2min(Idx2_inter)-S.chi2_bf,sin2T4(Idx2_inter),1,'spline')-S.sin2T4_bf;
%%
close all;
GetFigure;
plot(sin2T4,ones(numel(mNu4Sq),1),'k:','LineWidth',2);
hold on;
p1 = plot(sin2T4,chi2min-S.chi2_bf,'LineWidth',2);
hold off
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('\\chi^2 - \\chi^2_{min}'));
PrettyFigureFormat('FontSize',22);
xlim([1e-03 0.5]);
%%
leg = legend(p1,sprintf('Best fit: |{\\itU}_{e4}|^2 = %.3f (-%.3f +%.3f) eV^2',S.sin2T4_bf,sin2T4_errlow,sin2T4_errup),...
    'Location','northwest');
PrettyLegendFormat(leg);
print(gcf,'./plots/ksn2_ue4_bf_err.png','-dpng','-r450');

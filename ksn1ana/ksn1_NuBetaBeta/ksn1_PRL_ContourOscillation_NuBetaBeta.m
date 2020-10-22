% KSN1 prl plot: contours in oscillation parameter space
% Lisa, May 2020

% work flow for SterileAnalysis clas:
% 1. set up a (Multi-)RunAnalysis object "R" with general settings (Runlist, FSD,...)
% 2. set up SterileAnalysis object "S" with "R" as input argument
% 3. from here you can load chi2 grids / plot contours, change some basic settings
%% 1. RunAnalysis object
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
    'DataType','Real',...
    'FSDFlag','SibilleFull',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','OFF',...
    'ISCSFlag','Edep',...
    'TwinBias_Q',18573.73,...
    'SysBudget',24,...
    'pullFlag',99,...
    'NonPoissonScaleFactor',1};

R = MultiRunAnalysis(RunAnaArg{:});
R.chi2 = 'chi2CMShape';
%% 2. SterileAnalysis class
SterileArg = {'RunAnaObj',R,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',40};

S = SterileAnalysis(SterileArg{:});

%% 3. contour plot in oscillation parameter space (you can switch on/off foreign contours, all on by default
S.LoadGridFile('CheckSmallerN','ON','CheckLargerN','ON'); % if CheckSmallerN also look for grid with more/less nGridSteps
S.InterpMode = 'spline';           % waring: if contour is closed, spline interp sometimes sensitive to artefacts! Switch to "lin" in this case
S.Interp1Grid('RecomputeFlag','ON');% interpolate chi2 map -> nicer appearance of all plots. some

Arg = {'SavePlot','OFF','BestFit','OFF','FinalSensitivity','ON','Style','PRL'};
pHandle = S.ContourPlotOsci(Arg{:});

%%
hold on;
mbbGerda = 0.2;
mbbLegend200 = 0.06;
mbbLegend1000 = 0.03;
Ue4sq = logspace(-3,0,100); Ue4sq = Ue4sq(Ue4sq<=0.5);
sin22theta = 4 .* Ue4sq .* ( 1 - Ue4sq);
Dm41sq_Gerda = ( mbbGerda ./ Ue4sq).^2;
Dm41sq_Legend200 = ( mbbLegend200 ./ Ue4sq).^2;
Dm41sq_Legend1000 = ( mbbLegend1000 ./ Ue4sq).^2;
gerda = plot(sin22theta,Dm41sq_Gerda,':','LineWidth',2,'Color',rgb('HotPink'));

% hold on
% l200 = plot(sin22theta,Dm41sq_Legend200,'LineWidth',2);
% l1000 = plot(sin22theta,Dm41sq_Legend1000,'LineWidth',2);

%%
% normal hirachy
Hierarchy = 'IH';
close all
nSamples = 5e3;
mbbexpLim = 0.2; % experimental limit (GERDA)
if strcmp(Hierarchy,'NH')
    mbb3nuMax = 0.005; % normal hierarchy
else
    mbb3nuMax = 0.05; % inverted hierarchy
end

% randomize: 
phase    = 2.*rand(nSamples,1)-1;             % draw random phase between +-1 
m41      = 5.*rand(nSamples,1);          % masses between 0 and 50 eV
U4eSq    = 0.5.*rand(nSamples,1);             % mixing square between 0 and 0.5
m41Ue4Sq = m41.*U4eSq;%2.*rand(nSamples,1);               % draw random m41*Ue4^2
mbbexp   = abs((1-U4eSq).*mbb3nuMax+phase.*m41Ue4Sq);      % calculate experimentally measuread mbetabeta

% plot
h1 =hexscatter(mbbexp,m41Ue4Sq,'res',250);
%h1 =hexscatter(m41Ue4Sq,mbbexp,'res',250);

set(gcf,'Colormap',flipud(colormap('winter')));
%h1 = scatter(m41Ue4Sq,mbbexp,'o','MarkerEdgeColor','none','MarkerFaceColor',rgb('DodgerBlue'),'MarkerFaceAlpha',0.25,...
%    'SizeData',15);
hold on
%p1 = plot(linspace(0,max(m41Ue4Sq),10),mbbexpLim.*ones(10,1),'-','LineWidth',2,'Color',rgb('Black'));
p1 = plot(mbbexpLim.*ones(10,1),linspace(0,max(m41Ue4Sq),10),'-','LineWidth',2,'Color',rgb('Black'));

xlabel(sprintf('{\\itm}^{exp}_{\\beta\\beta}'))
ylabel(sprintf('|{\\itU}_{e4}|^{ 2} \\cdot {\\itm}_{41}'));
PrettyFigureFormat('FontSize',18);
x = linspace(0,max(m41Ue4Sq),10);
%pMax = plot(x,abs(mbb3nuMax+(1).*x),':','LineWidth',2.5,'Color',rgb('DodgerBlue'));
pMax = plot(abs(mbb3nuMax+(1).*x),x,':','LineWidth',2.5,'Color',rgb('DodgerBlue'));

%plot(m41Ue4Sq,abs(mbb3nuNH+(0.5).*m41Ue4Sq),':','LineWidth',2);
%plot(m41Ue4Sq,abs(mbb3nuNH+(0.25).*m41Ue4Sq),':','LineWidth',2);
hold off;

leg = legend([p1,h1,pMax], sprintf('{\\itm}_{\\beta\\beta}^{exp} < 0.2 eV'),...
    sprintf('Randomly drawn \\delta \\in [-1,1]'),...
    sprintf('\\delta = +1'),'Location','southeast');
legend boxoff;
leg.FontSize = 15;
t = title([sprintf('{\\itm}^{exp}_{\\beta\\beta} = | (1 - |{\\itU}_{e4}|^{ 2}) \\cdot {\\itm}_{\\beta\\beta}^{max}(3\\nu) + \\delta \\cdot |{\\itU}_{e4}|^{ 2} \\cdot {\\itm}_{41} | \n'),...
    sprintf('{\\itm}_{\\beta\\beta}^{max}(3\\nu) = %.1g eV (%s , non-degenerate)',mbb3nuMax,Hierarchy)],...
    'FontWeight','normal','FontSize',14);
xlim([0 1.2]);%max(m41Ue4Sq)]);
ylim([0 max(mbbexp(m41Ue4Sq<1.2))]);
%print(gcf,sprintf('./plots/NuBetaBeta%s.png',Hierarchy),'-dpng','-r300');

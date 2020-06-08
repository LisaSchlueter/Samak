% Lisa,  June 2020
% plot with neutrino mass distribution
% nuissance parameter free
% nuissance parameter + pull
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','mNu E0 Norm Bkg',...
    'DataType',DataType,...
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

T = MultiRunAnalysis(RunAnaArg{:});
T.chi2 = 'chi2CMShape';
%% settings sterile class
SterileArg = {'RunAnaObj',T,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',40};

S = SterileAnalysis(SterileArg{:});
%%
S.RunAnaObj.pullFlag = 13;
f = S.LoadGridFile;
mNuSq_2d = cellfun(@(x) x.par(1),f.FitResults);
E0_2d    = cellfun(@(x) x.par(2),f.FitResults);

mNuSq = reshape(cellfun(@(x) x.par(1),f.FitResults),[S.nGridSteps^2,1]);
E0    = reshape(cellfun(@(x) x.par(2),f.FitResults),[S.nGridSteps^2,1]);

%% 2D histogram: mNu
GetFigure
zlimMax = 15;
mNuSq_2d(abs(mNuSq_2d)>zlimMax)=NaN;
surf(S.sin2T4,S.mNu4Sq,mNuSq_2d,'EdgeColor','interp','FaceColor','interp');
PrettyFigureFormat('FontSize',24);
view([0 0 1])
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off
c = colorbar;
colormap('cool')
c.Label.String = sprintf('{\\itm}_\\nu^2 (eV^2)');
c.Label.FontSize = get(gca,'FontSize');
c.Limits=[-zlimMax zlimMax];
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
ylim([1 2e3])
xlim([1e-03 0.5])

if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag == 13
     title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free -  \\sigma({\\itE}_0) = 1 eV',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2); 
end

filename = sprintf('%smNuDistGrid_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');

%% 2D histogram: E0
GetFigure
surf(S.sin2T4,S.mNu4Sq,E0_2d+S.RunAnaObj.ModelObj.Q_i-18574,'EdgeColor','interp','FaceColor','interp');
PrettyFigureFormat('FontSize',24);
view([0 0 1])
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off
c = colorbar;
colormap('cool')
c.Label.String = sprintf('{\\itE}_0 - 18574 (eV)');
c.Label.FontSize = get(gca,'FontSize');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
ylim([1 2e3])
xlim([1e-03 0.5])

if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag == 13
     title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free -  \\sigma({\\itE}_0) = 1 eV',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2); 
end

filename = sprintf('%sE0DistGrid_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');


%% 1D histograms:
%% mnu
myFontSize = 24;
GetFigure;
h1 = histogram(mNuSq,'BinWidth',0.5);
PrettyFigureFormat('FontSize',myFontSize);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel('Occurrence');
xlim([-15 15])
if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag == 13
     title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free -  \\sigma({\\itE}_0) = 1 eV',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2); 
end

savedir = [getenv('SamakPath'),...
    sprintf('SterileAnalysis/plots/%s/%s/',...
    S.RunAnaObj.DataSet,S.RunAnaObj.DataType)];
filename = sprintf('%smNuDist_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');
%% E0
myFontSize = 24;
GetFigure;
h1 = histogram(E0+S.RunAnaObj.ModelObj.Q_i-18573.7,'BinWidth',0.05);
PrettyFigureFormat('FontSize',myFontSize);
xlabel(sprintf('{\\itE}_0 - 18573.7 (eV)'));
ylabel('Occurrence');
xlim([-0.5 3]);
if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
end

savedir = [getenv('SamakPath'),...
    sprintf('SterileAnalysis/plots/%s/%s/',...
    S.RunAnaObj.DataSet,S.RunAnaObj.DataType)];
filename = sprintf('%sE0Dist_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');
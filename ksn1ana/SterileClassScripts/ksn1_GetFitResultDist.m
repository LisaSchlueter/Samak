% Lisa,  June 2020
% plot with neutrino mass distribution
% nuissance parameter free
% nuissance parameter + pull
%% settings for runanalysis
DataType = 'Real';
%%
RunAnaArg = {'RunList','KNM1',...
    'fixPar','E0 Norm Bkg',...
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
S.RunAnaObj.pullFlag = 99;
f = S.LoadGridFile;
mNuSq_2d = cellfun(@(x) x.par(1),f.FitResults)';
E0_2d    = cellfun(@(x) x.par(2),f.FitResults)';

mNuSq = reshape(cellfun(@(x) x.par(1),f.FitResults),[S.nGridSteps^2,1]);
E0    = reshape(cellfun(@(x) x.par(2),f.FitResults),[S.nGridSteps^2,1]);

savedir = [getenv('SamakPath'),...
    sprintf('SterileAnalysis/plots/%s/%s/',...
    S.RunAnaObj.DataSet,S.RunAnaObj.DataType)];
%% 2D histogram: mNu
GetFigure
mNuSq_2d_plot = mNuSq_2d;
mNuSq_2d_plot(mNuSq_2d_plot<-5)=-5.5;
mNuSq_2d_plot(mNuSq_2d_plot>10)=11;

[C,h] = contourf(S.sin2T4,S.mNu4Sq,mNuSq_2d_plot,[-10 -2 -1 -0.5 0 1 3 10],...
    'Color',rgb('Black'),'LineWidth',1,'ShowText','ON');
clabel(C,h,'Color',rgb('Black'),'FontSize',12,...
    'FontSmoothing','on','FontWeight','bold','LabelSpacing',120,'Margin',1)
PrettyFigureFormat('FontSize',24);
view([0 0 1])
set(gca,'XScale','log');
set(gca,'YScale','log');
ylim([1 1e3])
xlim([1e-03 0.5])
c = colorbar;
colormap('parula')
c.Label.String = sprintf('{\\itm}_\\nu^2 (eV^2)');
c.Label.FontSize = get(gca,'FontSize');

xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
h.LevelList = [-6 -5 -2  -1.5 -1.1 -1  -0.5  0  1 5];
%c.Limits = [-10];
if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free (unconstrained)',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag == 13
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free -  \\sigma({\\itE}_0) = 1 eV',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif  S.RunAnaObj.pullFlag == 15
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif  S.RunAnaObj.pullFlag == 16
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 2 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif  S.RunAnaObj.pullFlag == 17
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 3 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
    
end

filename = sprintf('%smNuDistGrid_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');

%% 2D histogram: E0
GetFigure
E0_2d(E0_2d<-0.1) = -0.11;
[C,h] =contour(S.sin2T4,S.mNu4Sq,E0_2d,[-10 -2 -1 -0.5 0 1 3 10],...
    'Color',rgb('Black'),'LineWidth',1.5,'ShowText','ON');
clabel(C,h,'Color',rgb('Black'),'FontSize',12,...
    'FontSmoothing','on','FontWeight','bold','LabelSpacing',180)
h.Fill = 'on';
%surf(S.sin2T4,S.mNu4Sq,E0_2d+S.RunAnaObj.ModelObj.Q_i-18574,'EdgeColor','interp','FaceColor','interp');
PrettyFigureFormat('FontSize',24);
%view([0 0 1])
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off
c = colorbar;
colormap('parula')
c.Label.String = sprintf('{\\itE}_0 - %.1f (eV)',T.ModelObj.Q_i);
c.Label.FontSize = get(gca,'FontSize');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
%ylim([1 2e3])
%xlim([1e-03 0.5])
ylim([1 1e2])
xlim([1e-03 0.5])
h.LevelList = [-0.2 -0.1,0,0.05,0.1,0.25,0.5,1];

if S.RunAnaObj.pullFlag==12
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1.94 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag==99
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free (unconstrained)',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif S.RunAnaObj.pullFlag == 13
    title(sprintf('%.0f eV range - {\\itm}_\\nu^2 free -  \\sigma({\\itE}_0) = 1 eV',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif  S.RunAnaObj.pullFlag == 15
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 1 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
elseif  S.RunAnaObj.pullFlag == 16
    title(sprintf('%.0f eV range - \\sigma({\\itm}_\\nu^2) = 2 eV^2',S.range),...
        'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
end

filename = sprintf('%sE0DistGrid_%.0feV_pull%.0f.png',...
    savedir,S.range,S.RunAnaObj.pullFlag);
print(gcf,filename,'-dpng','-r400');



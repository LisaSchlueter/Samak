% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + pull
% fixed parameter
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
    'pullFlag',12,...
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

S = St% Test SterileAnalysis class
% Lisa, May2020
% plot with neutrino mass as
% nuissance parameter free
% nuissance parameter + pull
% fixed parameter
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
f = S.LoadGridFile;
mNuSq = reshape(cellfun(@(x) x.par(1),f.FitResults),[50*50,1]);
E0    = reshape(cellfun(@(x) x.par(2),f.FitResults),[50*50,1]);


%% mnu
myFontSize = 24;
GetFigure;
h1 = histogram(mNuSq);
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
h1 = histogram(E0+S.RunAnaObj.ModelObj.Q_i-18573.7);
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
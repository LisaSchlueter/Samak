% plot script for KSN1 sterile contour plot
% Lisa, April 2020
%% settings
%Fitrium = 'ON';
CL = 95;%, 95 99];
SavePlot = 'ON';
PlotSplines = 'ON'; % plot smooth spline instead of points -> fails for closed contours
range = [95:-5:45,41,40];%[65:-5:45,41,40];%[95,90,80,75,70];%
nGridSteps = 50;
chi2Str = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
SysBudget = 24;
%% load grid (or calculate if doesn't exist)
nContours = numel(range);
mnu4Sq   = cell(nContours,1);
sin2T4   = cell(nContours,1);
chi2     = cell(nContours,1);
chi2_ref = cell(nContours,1);

for i=1:numel(range)
    GridArg = {'range',range(i),...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
     'SysBudget',SysBudget};
[mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile] = KSN1GridSearch(...
    'DataType',DataType,GridArg{:});
end
%% collect plot argument

if strcmp(chi2Str,'chi2Stat')
    chi2Label = 'stat. only';
else
    chi2Label = 'stat. and syst.';
end
titleStr = sprintf('%s , %s at %.0f%% C.L.',DataType,chi2Label,CL);

LineStyles = {'-','-.',':','-','--','-.',':','--','-','-',':','-.','-'};
Colors = {'DodgerBlue','Orange','DarkSlateGray','FireBrick','Magenta','LimeGreen','CadetBlue',...
    'Navy','ForestGreen','PowderBlue','Pink','DarkOrange','Black'};
legStr = cell(nContours,1);
pl = cell(nContours,1);
for i=1:nContours
    PlotArg ={'mnu4Sq',mnu4Sq{i},...
        'sin2T4',sin2T4{i},...
        'chi2',chi2{i},'chi2_ref',chi2_ref{i},...
        'CL',CL,...
        'titleStr',titleStr,'LineStyle',LineStyles{i},...
        'Color',Colors{i}};
    %% plot
    if strcmp(DataType,'Twin')
        Method = 'Old';
    else
        Method = 'New';
    end
if i>1
    [pl{i},~] = KSN1ContourPlot(PlotArg{:},'PlotSplines',PlotSplines,'HoldOn','ON','Method',Method);
else
    [pl{i},~] = KSN1ContourPlot(PlotArg{:},'PlotSplines',PlotSplines,'Method',Method);
end

legStr{i} = sprintf('%.0f eV range',range(i));
end

leg = legend([pl{:}]',legStr{:},'EdgeColor',rgb('Silver'),'Location','southwest');
leg.Title.String = 'Lower fit boundary';
leg.Title.FontWeight = 'normal';
if numel(range)>5
    leg.NumColumns=2;
end
ylim([1 1e4]);
if strcmp(DataType,'Real')
    xlim([1e-03 0.5]);
else
    xlim([4e-03 0.5]);
end
plotdir = [getenv('SamakPath'),'ksn1ana/LisaSterile/plots/'];
plotname = sprintf('%sksn1qUScan_%s.png',plotdir,DataType);
print(gcf,plotname,'-dpng','-r450');




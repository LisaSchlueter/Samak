% investigate impact of different grid types on contour 
nGridSteps = 30;
mNu4SqTestGrid = [4,5,3];
Mode = 'Compute';
%% configure RunAnalysis object
chi2 = 'chi2Stat';
DataType = 'Twin';
range = 40;
InterpMode = 'spline';
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
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
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
%%
%S.GridSearch('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
InterpMode = 'spline';
S.InterpMode = InterpMode;
switch Mode
    case 'Compute'
        for i=1:numel(mNu4SqTestGrid)
            S.GridSearch('mNu4SqTestGrid',mNu4SqTestGrid(i),'ExtmNu4Sq','ON');
        end
        S.GridSearch('mNu4SqTestGrid','OFF','ExtmNu4Sq','ON');
    case 'Display'
        Color = {'DodgerBlue','Orange','FireBrick','SkyBlue','LimeGreen'};
        LineStyle = {'--',':','-.','-.','-.',':','--','-'};
        pHandle = cell(numel(mNu4SqTestGrid)+4,1);
        legStr = cell(numel(mNu4SqTestGrid)+4,1);
        
        % 50 grid
        S.nGridSteps= 50;  
        S.LoadGridFile('mNu4SqTestGrid','OFF','ExtmNu4Sq','ON','IgnoreKnm2FSDbinning','ON','CheckLargerN','OFF');
        mNu50 = S.mNu4Sq(:,1);
        S.Interp1Grid('Maxm4Sq',38.2^2)
        pHandle{1} = S.ContourPlot('HoldOn','OFF','Color',rgb('Black'));
        legStr{1} = sprintf('50 \\times 50 log-log grid (%.0f points for {\\itm}_4^2\\geq900 eV^2)',sum(mNu50>=29.999^2));
       
        % 30 grid
        S.nGridSteps= 30;
        S.LoadGridFile('mNu4SqTestGrid','OFF','CheckExtmNu4Sq','ON','CheckLargerN','OFF','IgnoreKnm2FSDbinning','ON','ExtmNu4Sq','ON');
        mNu30 = S.mNu4Sq(:,1);
        S.Interp1Grid('Maxm4Sq',38^2)
        pHandle{2} = S.ContourPlot('HoldOn','ON','Color',rgb('Black'),'LineStyle','-.');
        legStr{2} = sprintf('30 \\times 30 log-log grid (%.0f points for {\\itm}_4^2\\geq900 eV^2)',sum(mNu30>=29.999^2));
     
        S.LoadGridFile('mNu4SqTestGrid',2,'ExtmNu4Sq','ON','IgnoreKnm2FSDbinning','ON','CheckLargerN','OFF','ExtmNu4Sq','ON');
        mNu30 = S.mNu4Sq(:,1);
        S.Interp1Grid('Maxm4Sq',38.2^2)
        pHandle{3} = S.ContourPlot('HoldOn','ON','Color',rgb(Color{1}),'LineStyle','-.');
        legStr{3} = sprintf('30 \\times 30 with %.0f points for {\\itm}_4^2\\geq900 eV^2',sum(mNu30>=29.999^2));
     
        % 25 grid
        S.nGridSteps= 25;   
        S.LoadGridFile('mNu4SqTestGrid','OFF','ExtmNu4Sq','OFF','IgnoreKnm2FSDbinning','ON','CheckLargerN','OFF','ExtmNu4Sq','ON');
        mNu25reg = S.mNu4Sq(:,1);
        S.Interp1Grid('Maxm4Sq',38.2^2)
        pHandle{4} = S.ContourPlot('HoldOn','ON','Color',rgb('Black'),'LineStyle',':');
        legStr{4} = sprintf('25 \\times 25 log-log grid (%.0f points for {\\itm}_4^2\\geq900 eV^2)',sum(mNu25reg>=(29.999)^2));
       
       
        for i=1:numel(mNu4SqTestGrid)
            S.LoadGridFile('mNu4SqTestGrid',mNu4SqTestGrid(i),'ExtmNu4Sq','ON','IgnoreKnm2FSDbinning','ON','CheckLargerN','OFF');
             mNu25 =S.mNu4Sq(:,1);
            S.Interp1Grid('Maxm4Sq',38.2^2)
            pHandle{i+4} = S.ContourPlot('HoldOn','ON','Color',rgb(Color{i+1}),'LineStyle',LineStyle{i});
            legStr{i+4} = sprintf('25 \\times 25 with %.0f points for {\\itm}_4^2\\geq900 eV^2',sum(mNu25>=29.999^2));
        end
end
%%
%leg = legend([pHandle{:}],legStr);
OrderIdx = [4,2,1,5,3,6];
leg = legend([pHandle{OrderIdx}],legStr{OrderIdx});
PrettyLegendFormat(leg);

%% 
plotdir = [getenv('SamakPath'),'ksn2ana/ksn2_GridSize/plots/'];
plotname = sprintf('%sksn2_DifferentGrids_%s.png',plotdir,InterpMode);
MakeDir(plotdir);
print(plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);
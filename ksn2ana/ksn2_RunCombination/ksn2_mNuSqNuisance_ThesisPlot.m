% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
% as in ksn2 paper!
%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';%'Real';
nGridSteps = 40;
range = 40;
Mode = 'Compute';
LocalFontSize = 19;
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','mNu E0 Norm Bkg',...%free par
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
    'range',range,...
    'LoadGridArg',{'mNu4SqTestGrid',5}};

S = SterileAnalysis(SterileArg{:});
S.InterpMode = 'spline';
%
S.LoadGridFile(S.LoadGridArg{:},'ExtmNu4Sq','ON'); 
S.Interp1Grid('nInter',1e3);
S.FindBestFit;
mNu4Sq_bf = S.mNu4Sq_bf;
sin2T4_bf = S.sin2T4_bf;

PlotPar = S.mNuSq;

if strcmp(DataType,'Real')
    ContourVec = [-5 -1,0,0.28,1,2,10];
    LabelSpacing = 380;
      BF = 'ON';
else
    ContourVec = [-5 -1,-0.3 0,0.3,1,2,10];
    LabelSpacing = 10000;
    BF = 'OFF';
end
Splines = 'ON';
%f = GetFigure;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
%f.Renderer = 'painters';
%
for i=1:numel(ContourVec)
    [M,p1] = contour3(S.sin2T4,S.mNu4Sq,PlotPar,[ContourVec(i),ContourVec(i)],...
        'Color',rgb('LightSalmon'),...
        'ShowText','off','LineWidth',2.5,...
          'LabelSpacing',LabelSpacing,...
          'LineStyle',':');
    hold on;
   % cl = clabel(M,p1,'FontSize',18,'FontName','Times New Roman');
end
%
[pFree,pFix] = S.PlotmNuSqOverview('PullmNuSq','OFF','SavePlot','OFF','HoldOn','ON','BestFit','OFF');
pFree.LineColor = rgb('FireBrick');
if strcmp(BF,'ON')
   pbf =  plot3(sin2T4_bf,mNu4Sq_bf,99,...
       'p','MarkerIndices',[1,1],'MarkerFaceColor','red',...
       'LineWidth',1,'Color',pFree.Color,'MarkerSize',17,'MarkerFaceColor',pFree.Color);
end
view(2)
grid off
pFix.LineWidth = 3;
pFix.LineStyle = '-.';
pFree.LineWidth = 3;

%
pNone = plot3(NaN,NaN,NaN,'w','LineWidth',0.1);
leg = legend([pNone,pFix,pFree,p1],...
    'KNM2',...
    sprintf('I)  Fixed {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('II) Free {\\itm}_\\nu^2 unconstrained'),...
    sprintf('Isoline: {\\itm}_\\nu^2 (eV^2) best fit for II)'),...
    'EdgeColor',rgb('Silver'),'Location','southwest','box','off');
if strcmp(DataType,'Twin')
    legend boxon
    PrettyLegendFormat(leg,'alpha',0.5);
end
%

title('');
%
%PRLFormat;set(gca,'FontSize',LocalFontSize);
PrettyFigureFormat('FontSize',LocalFontSize);
leg.FontSize = LocalFontSize-2;
leg.Position(2) = 0.24;
ax = gca;
% ax.XLabel.FontSize = LocalFontSize;
% ax.YLabel.FontSize = LocalFontSize;

ylim([0.1 2e3]);
yticks([0.1 1 10 100 1e3])
xlim([4e-03 0.5]);
set(gca,'YScale','log')
set(gca,'XScale','log')
% make nice labels by hand
if strcmp(DataType,'Real')
    BoxArg = { 'Color','Black','EdgeColor','none','FontSize',LocalFontSize-2,'FontName',get(gca,'FontName')};
    if exist('t1','var')
    t1.delete;
    tx1.delete;
    t2.delete;
    tx2.delete;
    t3.delete;
    tx3.delete;
    t4.delete;
    tx4.delete;
    t5.delete;
    tx5.delete;
    t6.delete;
    tx6.delete;
    t7.delete;
    tx7.delete;
    t8.delete;
    tx8.delete;
        t9.delete;
    tx9.delete;
    end
    t1 = annotation('rectangle',[0.82,0.2,0.07,0.05],'FaceColor',rgb('White'),'Color','none');
    tx1 = annotation('textbox',[0.82,0.27,0.05,0.001],'String',0.28,BoxArg{:});
    
    t2 = annotation('rectangle',[0.84,0.3,0.03,0.04],'FaceColor',rgb('White'),'Color','none');
    tx2 = annotation('textbox',[0.84,0.35,0.05,0.005],'String',0,BoxArg{:});

    t3 = annotation('ellipse',[0.83,0.37,0.04,0.05],'FaceColor',rgb('White'),'Color','none'); 
    tx3 = annotation('textbox',[0.835,0.42,0.01,0.01],'String','-1',BoxArg{:});
 
    t4 = annotation('ellipse',[0.8,0.57,0.047,0.06],'FaceColor',rgb('White'),'Color','none');
    tx4 = annotation('textbox',[0.8,0.585,0.05,0.05],'String','-5',BoxArg{:});
  
    t5 = annotation('rectangle',[0.732,0.7,0.04,0.04],'FaceColor',rgb('White'),'Color','none');
    tx5 = annotation('textbox',[0.73,0.72,0.05,0.03],'String','10',BoxArg{:});
  
    t6 = annotation('rectangle',[0.53,0.81,0.03,0.04],'FaceColor',rgb('White'),'Color','none');
    tx6 = annotation('textbox',[0.53,0.82,0.05,0.05],'String','2',BoxArg{:});
  
    t7 = annotation('rectangle',[0.425,0.82,0.02,0.03],'FaceColor',rgb('White'),'Color','none');
    tx7 = annotation('textbox',[0.42,0.82,0.05,0.05],'String','1',BoxArg{:});   
 
    t8 = annotation('rectangle',[0.302,0.62,0.07,0.05],'FaceColor',rgb('White'),'Color','none');
    tx8 = annotation('textbox',[0.305,0.64,0.05,0.02],'String','0.28',BoxArg{:});
   
    t9 = annotation('rectangle',[0.2,0.88,0.07,0.04],'FaceColor',rgb('White'),'Color','none');
    tx9 = annotation('textbox',[0.2,0.874,0.05,0.05],'String','0.28',BoxArg{:});
end

%% save
axes.SortMethod='ChildOrder';
name_i = strrep(S.DefPlotName,'_mNuE0BkgNorm','');
plotname = sprintf('%s_mNuSqOverviewmNuSq_%.2gCL.png',name_i,S.ConfLevel);
%print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
%
%ylabel(sprintf('{\\itm}_4^2 (eV^{2})'));
export_fig(strrep(plotname,'.png','.pdf'),'-opengl');
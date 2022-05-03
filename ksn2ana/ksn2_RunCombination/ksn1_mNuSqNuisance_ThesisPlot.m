% ksn2 calculate chi2 grid search
% compare m2 free, m2 nuisance parameter
% as in ksn2 paper!
%% settings that might change
DataType              = 'Real';
range                 = 40;
chi2                  = 'chi2CMShape';
freePar               = 'mNu E0 Norm Bkg';
LocalFontSize = 19;
%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.064;
end
Real = MultiRunAnalysis('RunList','KNM1',...
    'chi2',chi2,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'SysBudget',24,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'AngularTFFlag','OFF',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','OFF',...
    'BKG_PtSlope',0);
Real.exclDataStart = Real.GetexclDataStart(range);
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

S = SterileAnalysis(SterileArg{:});


Real.AngularTFFlag = 'ON';
Real.FSDFlag = 'KNM2_0p1eV';
Real.ModelObj.BKG_PtSlope = -2.2*1e-06;
Real.ELossFlag = 'KatrinT2A20';
S.RunAnaObj.SysBudget = 200;
S.nGridSteps = 50;


S.LoadGridArg = {'mNu4SqTestGrid',2};
S.InterpMode = 'Mix';
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('Maxm4Sq',36^2);
S.FindBestFit;
mNu4Sq_bf = S.mNu4Sq_bf;
sin2T4_bf = S.sin2T4_bf;

PlotPar = S.mNuSq;

if strcmp(DataType,'Real')
    ContourVec = [-10 -5.3 -2 -1,0,1,2,10];
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
%
[pFree,pFix] = S.PlotmNuSqOverview('PullmNuSq','OFF','SavePlot','OFF','HoldOn','ON','BestFit','OFF');
pFree.LineColor = rgb('FireBrick');
if strcmp(BF,'ON')
   pbf =  plot3(sin2T4_bf,mNu4Sq_bf,99,...
       'p','MarkerIndices',[1,1],'MarkerFaceColor','red',...
       'LineWidth',1,'Color',pFree.Color,'MarkerSize',17,'MarkerFaceColor',pFree.Color);
   
   S.FindBestFit;
    pbffix =  plot3(S.sin2T4_bf,S.mNu4Sq_bf,99,...
       'h','MarkerIndices',[1,1],'MarkerFaceColor',rgb('DodgerBlue'),...
       'LineWidth',1,'Color',pFix.Color,'MarkerSize',17,'MarkerFaceColor',pFix.Color);
end
view(2)
grid off
pFix.LineWidth = 3;
pFix.LineStyle = '-.';
pFree.LineWidth = 3;

%
pNone = plot3(NaN,NaN,NaN,'w','LineWidth',0.1);
leg = legend([pNone,pFix,pFree,p1],...
    'KNM1',...
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
leg.Position(2) = 0.2;

%leg.Title.String = 'KNM1';
%leg.Title.FontWeight = 'normal';

ylim([1 1600]);
yticks([0.1 1 10 100 1e3])
xlim([8e-03 0.5]);
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
    t1 = annotation('rectangle',[0.52,0.335,0.05,0.05],'FaceColor',rgb('White'),'Color','none');
    tx1 = annotation('textbox',[0.535,0.4,0.05,0.001],'String',-2,BoxArg{:});
    
    t2 = annotation('rectangle',[0.72,0.395,0.05,0.04],'FaceColor',rgb('White'),'Color','none');
    tx2 = annotation('textbox',[0.71,0.44,0.05,0.005],'String',-5.3,BoxArg{:});

    t3 = annotation('ellipse',[0.81,0.39,0.05,0.05],'FaceColor',rgb('White'),'Color','none'); 
    tx3 = annotation('textbox',[0.805,0.434,0.01,0.01],'String',-10,BoxArg{:});

    t5 = annotation('rectangle',[0.605,0.7,0.03,0.04],'FaceColor',rgb('White'),'Color','none');
    tx5 = annotation('textbox',[0.6,0.72,0.05,0.03],'String','10',BoxArg{:});
  
    t6 = annotation('rectangle',[0.37,0.75,0.03,0.04],'FaceColor',rgb('White'),'Color','none');
    tx6 = annotation('textbox',[0.375,0.75,0.05,0.05],'String','2',BoxArg{:});
  
    t7 = annotation('rectangle',[0.354,0.79,0.02,0.03],'FaceColor',rgb('White'),'Color','none');
    tx7 = annotation('textbox',[0.35,0.79,0.05,0.05],'String','1',BoxArg{:});   
 
    t8 = annotation('rectangle',[0.25,0.58,0.035,0.05],'FaceColor',rgb('White'),'Color','none');
    tx8 = annotation('textbox',[0.25,0.615,0.05,0.02],'String',-1,BoxArg{:});
   
    t9 = annotation('rectangle',[0.2,0.86,0.035,0.04],'FaceColor',rgb('White'),'Color','none');
    tx9 = annotation('textbox',[0.2,0.86,0.05,0.05],'String',-1,BoxArg{:});
    
    t4 = annotation('rectangle',[0.23,0.68,0.035,0.04],'FaceColor',rgb('White'),'Color','none');
    tx4 = annotation('textbox',[0.235,0.68,0.05,0.05],'String',0,BoxArg{:});
end

% save
axes.SortMethod='ChildOrder';
name_i = strrep(S.DefPlotName,'_mNuE0BkgNorm','');
plotname = sprintf('%s_mNuSqOverviewmNuSq_%.2gCL.png',name_i,S.ConfLevel);
%print(gcf,plotname,'-dpng','-r450');
fprintf('save plot to %s \n',plotname);
%
%ylabel(sprintf('{\\itm}_4^2 (eV^{2})'));
export_fig(strrep(plotname,'.png','.pdf'),'-opengl');
%% compare KSN-2 numass with KNM-2
mNuSq_knm2    = 0.26; % central values from publication
mNuSqErr_knm2 = 0.34; % 1 sigma from publication main result
DataType = 'Real';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_mNuSqPulls/results/'];
MakeDir(savedir);
savefile = sprintf('%sksn2_ComparemNuSq_KNM2_%s_mNuSq%.2feV2_mNuSqErr%-.2feV2.mat',savedir,DataType,mNuSq_knm2,mNuSqErr_knm2);
if exist(savefile,'file')
    load(savefile);
    fprintf('load from file %s \n',savefile)
else
    
    %% settings that might change
    chi2 = 'chi2CMShape';
    nGridSteps = 40;
    range = 40;
    Mode = 'Compute';
    % configure RunAnalysis object
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
        'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};
    
    S = SterileAnalysis(SterileArg{:});
    S.InterpMode = 'spline';
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('nInter',1e3);
    S.ContourPlot('BestFit','ON'); close;
    
    sin2T4         = S.sin2T4;
    mNu4Sq         = S.mNu4Sq;
    mNuSq          =  S.mNuSq;
    sin2T4_contour = S.sin2T4_contour;
    mNu4Sq_contour = S.mNu4Sq_contour;
    sin2T4_bf      = S.sin2T4_bf;
    mNu4Sq_bf      = S.mNu4Sq_bf;
    
    ContourVec = [mNuSq_knm2-mNuSqErr_knm2,mNuSq_knm2,mNuSq_knm2+mNuSqErr_knm2];
    
    GetFigure
    [M1,p1] = contour3(sin2T4,mNu4Sq,mNuSq,[ContourVec(1),ContourVec(1)]); hold on;%close;
    sin2T4_iso1 = M1(1,2:M1(2,1));
    mNu4Sq_iso1 = M1(2,2:M1(2,1));
    [M3,p1] = contour3(sin2T4,mNu4Sq,mNuSq,[ContourVec(3),ContourVec(3)]); %close;
    sin2T4_iso3 = M3(1,2:M3(2,1));
    mNu4Sq_iso3 = M3(2,2:M3(2,1));
    
    save(savefile,'sin2T4','mNu4Sq','sin2T4_contour','mNu4Sq_contour',...
        'sin2T4_bf','mNu4Sq_bf','sin2T4_iso1','sin2T4_iso3','ContourVec','mNu4Sq_iso1','mNu4Sq_iso3');
    fprintf('save to file %s \n',savefile)
end

if strcmp(DataType,'Real')
    DataStr = 'exclusion';
else
    DataStr = 'sensitivity';
end


%%
GetFigure;
%plot(sin2T4_iso1,mNu4Sq_iso1,'LineWidth',2)
%hold on;
%plot(sin2T4_iso3,mNu4Sq_iso3,'LineWidth',2)
[l1,a1]=boundedline(sin2T4_iso1,mNu4Sq_iso1,[sin2T4_iso1-1e-05; 1e-05.*ones(1,numel(sin2T4_iso1))]','orientation','horiz');
hold on;
[l1tmp,a1tmp] = boundedline(0.5.*ones(10,1),linspace(max(mNu4Sq_iso1),1600,10),[0.5-1e-05;0],'orientation','horiz');
l1tmp.delete;
[l2tmp,a2tmp] = boundedline(0.5.*ones(10,1),linspace(0.1,min(mNu4Sq_iso1),10),[0.5-1e-05;0],'orientation','horiz');
l2tmp.delete;
[l3,a3]=boundedline(sin2T4_iso3,mNu4Sq_iso3,[1e-05,1],'orientation','horiz');
a3.FaceColor = rgb('White');
pContour = plot(sin2T4_contour,mNu4Sq_contour,'-k','LineWidth',3,'Color',rgb('ForestGreen'));
if strcmp(DataType,'Real')
    pbf = plot(sin2T4_bf,mNu4Sq_bf,'x','LineWidth',2,'MarkerSize',7,'Color',pContour.Color);
end
% style
a1.FaceColor = rgb('LightGray');l1.Color = rgb('SlateGray');
a1tmp.FaceColor = rgb('LightGray');l3.Color = rgb('SLateGray');
a2tmp.FaceColor = rgb('LightGray');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([1e-03 0.5;])
ylim([0.1 1600])
PrettyFigureFormat;
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU}_{e4}|^2'));
leg = legend([pContour,a1tmp],sprintf('KATRIN %s at 95%% C.L.: {\\itm}_\\nu^2 free',DataStr),...
    sprintf('{\\itm}_\\nu^2 \\in %.2f \\pm %.2f eV^2 (KNM-2)',mNuSq_knm2,mNuSqErr_knm2),...
    'Location','southwest');
PrettyLegendFormat(leg);

%% save plot
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_mNuSqPulls/plots/'];
MakeDir(pltdir);
pltname = strrep(strrep(savefile,'results','plots'),'.mat','.png');
print(gcf,pltname,'-dpng');
  fprintf('save plot  to  %s \n',pltname)
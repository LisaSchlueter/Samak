% perform grid search in vicinity of best fit obtained from coarse fit

%% settings that might change
chi2 = 'chi2CMShape';
DataType = 'Real';
nGridSteps_FullGrid   = 50;
nGridSteps_BfImp = 40;

freePar = 'mNuSq E0 Norm Bkg';

savedir = sprintf('%sksn2ana/ksn2_BestFit/results/',getenv('SamakPath'));
savefile = sprintf('%sksn2_ImpBestFit_%s_%s_%s_nGridStepsFull%.0f_nGridStepsImp%.0f.mat',...
    savedir,DataType,strrep(freePar,' ',''),chi2,nGridSteps_FullGrid,nGridSteps_BfImp);
if exist(savefile,'file')
    load(savefile)
else  
    LoadGridArg = {'ExtmNu4Sq','OFF','mNu4SqTestGrid',5};
    range = 40;
    %% configure RunAnalysis object
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
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps_FullGrid,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',LoadGridArg};
    
    S = SterileAnalysis(SterileArg{:});
    %% load full grid
    S.LoadGridFile(S.LoadGridArg{:});
    
    % find grid point with smallest chi^2 (no interpolation!)
    [row, col]    = find(S.chi2 == min(S.chi2(:)));
    mNu4SqFull_bf = S.mNu4Sq(col,1);
    sin2T4Full_bf = S.sin2T4(1,row);
    chi2Full_bf = min(S.chi2(:));
    
    % construct grid around best fit
    S.nGridSteps = nGridSteps_BfImp;
    mNu4Sq_ImpGrid = logspace(log10(S.mNu4Sq(col-1,1)),log10(S.mNu4Sq(col+1,1)),S.nGridSteps)'; % has to be nGridSteps x 1
    sin2T4_ImpGrid = logspace(log10(S.sin2T4(1,row-1)),log10(S.sin2T4(1,row+1)),S.nGridSteps); % has to be 1 x nGridSteps
    ExternalGrid = {mNu4Sq_ImpGrid,sin2T4_ImpGrid};
    
    % do interpolate for comparison
    S.Interp1Grid;
    S.FindBestFit;
    mNu4SqFull_bf_inter = S.mNu4Sq_bf;
    sin2T4Full_bf_inter = S.sin2T4_bf;
    chi2Full_bf_inter = S.chi2_bf;
    
    %% load grid in vicinity of best fit
    S.nGridSteps = nGridSteps_BfImp;
 %   S.GridSearch('ExtGrid',ExternalGrid);
    S.LoadGridFile('ExtGrid',ExternalGrid);
    [row, col]    = find(S.chi2 == min(S.chi2(:)));
    mNu4Sq_bf = S.mNu4Sq(col,1);
    sin2T4_bf = S.sin2T4(1,row);
    chi2_bf = min(S.chi2(:));
    
    S.Interp1Grid('Minm4Sq',min(mNu4Sq_ImpGrid),'Maxm4Sq',max(mNu4Sq_ImpGrid));
    S.FindBestFit;
    mNu4Sq_bf_inter = S.mNu4Sq_bf;
    sin2T4_bf_inter = S.sin2T4_bf;
    chi2_bf_inter   = S.chi2_bf;
    mNu4SqGrid_inter = S.mNu4Sq;
    sin2T4Grid_inter = S.sin2T4;
    chi2Grid_inter = S.chi2;
    mNuSq_bf = S.mNuSq_bf;
    E0_bf = S.E0_bf;
    
    [DeltaChi2, SignificanceBF] = S.CompareBestFitNull;
    SigmaBF = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF,'nPar',2);
    chi2_Null = S.chi2_Null;
    dof = S.dof;
    pVal = 1-chi2cdf(chi2_bf_inter,dof);
    %%
    DefPltTitle = S.GetPlotTitle;
    MakeDir(savedir);
    save(savefile,...
        'ExternalGrid','sin2T4Grid_inter','mNu4SqGrid_inter','chi2Grid_inter',...%small grid
        'mNu4Sq_bf_inter','sin2T4_bf_inter','chi2_bf_inter',...                  %small grid
        'mNu4Sq_bf','sin2T4_bf','chi2_bf',...                                    %small grid
        'mNu4SqFull_bf','sin2T4Full_bf','chi2Full_bf',...                        %full large grid
        'mNu4SqFull_bf_inter','sin2T4Full_bf_inter','chi2Full_bf_inter',...      %full large grid
        'DefPltTitle','RunAnaArg','SterileArg',...
        'DeltaChi2', 'SignificanceBF','SigmaBF','chi2_Null','pVal','dof','E0_bf','mNuSq_bf');
end

GetFigure;
surf(sin2T4Grid_inter,mNu4SqGrid_inter,chi2Grid_inter,'EdgeColor','interp','FaceColor','interp');
view(2);
set(gca,'XScale','log');
set(gca,'YScale','log');
c = colorbar;
hold on;
p1 = plot3(sin2T4Full_bf,mNu4SqFull_bf,99,'x','MarkerSize',10,'LineWidth',2,'Color',rgb('White'));
p1_inter = plot3(sin2T4Full_bf_inter,mNu4SqFull_bf_inter,99,'o','MarkerSize',10,'LineWidth',2,'Color',rgb('White'));
p2 = plot3(sin2T4_bf,mNu4Sq_bf,99,'x','MarkerSize',10,'LineWidth',2,'Color',rgb('IndianRed'));
p2_inter = plot3(sin2T4_bf_inter,mNu4Sq_bf_inter,99,'o','MarkerSize',10,'LineWidth',2,'Color',rgb('IndianRed'));

xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat;

leg = legend([p1,p1_inter,p2,p2_inter],...
    sprintf('Best fit grid point ({\\itm}_4^2 \\in [1,1600] eV^2 and |{\\itU}_{e4}|^2 \\in [10^{-3},0.5])'),...
    sprintf('2D spline interp.  ({\\itm}_4^2 \\in [1,1600] eV^2 and |{\\itU}_{e4}|^2 \\in [10^{-3},0.5])'),...
    sprintf('Best fit grid point ({\\itm}_4^2 \\in [%.0f,%.0f] eV^2 and |{\\itU}_{e4}|^2 \\in [%.3f,%.3f])',min(ExternalGrid{1}),max(ExternalGrid{1}),min(ExternalGrid{2}),max(ExternalGrid{2})),...
    sprintf('2D spline interp.  ({\\itm}_4^2 \\in [%.0f,%.0f] eV^2 and |{\\itU}_{e4}|^2 \\in [%.3f,%.3f])',min(ExternalGrid{1}),max(ExternalGrid{1}),min(ExternalGrid{2}),max(ExternalGrid{2})),'Location','southwest');

title(DefPltTitle,'FontWeight','normal','FontSize',get(gca,'FontSize'));
PrettyLegendFormat(leg);
leg.Title.String = 'Location of best fit';
leg.Title.FontWeight = 'normal';
grid off
c.Label.String = sprintf('\\chi^2');
ax = gca;
c.Label.FontSize = ax.XLabel.FontSize;

pltdir = strrep(savedir,'results','plots');
% save
MakeDir(pltdir);
pltname = strrep(strrep(savefile,'results','plots'),'.mat','.png');
print(gcf,pltname,'-dpng','-r300')
fprintf('save plot to %s \n',pltname);
%%
%
fprintf('Best fit from fine grid: \n');
fprintf('sint4^2 = %.4f, m4^2 = %.2f eV^2 chi2min = %.2f, p = %.2f (%.0f dof)\n',sin2T4_bf_inter,mNu4Sq_bf_inter,chi2_bf_inter,pVal,dof);
fprintf('fit improves over null hypothesis by = %.2f  %.1f%% C.L. or %.2f sigma \n', DeltaChi2,SignificanceBF*100,SigmaBF);
fprintf('mNuSq = %.2f eV^2 ,E0 = %.1f eV \n', mNuSq_bf,E0_bf);


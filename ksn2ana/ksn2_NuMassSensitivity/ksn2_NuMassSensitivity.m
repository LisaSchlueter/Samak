% impact of 3+1 model extionsion on neutrino mass sensitivity
% no external constraints
chi2 = 'chi2CMShape';
DataType = 'Twin';%'Real';%'Twin';
SavePlt = 'ON';
%mNuSq = -2:0.05:1;
if strcmp(DataType,'Real')
mNuSq = -1:0.01:2.5;
else
mNuSq = -1:0.01:2.5;
end
%mNuSq = -2:0.1:1;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_NuMassSensitivity_%s_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.mat',...
    savedir,DataType,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));

if exist(savefile,'file') && 1==2
    load(savefile,'mNuSq','chi2min');
else
    %% settings that might change
    range = 40;
    if strcmp(DataType,'Real')
        nGridSteps = 40;%50;
        LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
    else
        nGridSteps = 40;
        LoadGridArg =   {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
    end
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
        'LoadGridArg',LoadGridArg};
    
    S = SterileAnalysis(SterileArg{:});
    S.LoadGridFile(S.LoadGridArg{:});
    S.InterpMode = 'spline';
    S.Interp1Grid('Maxm4Sq',38^2);
    
    S.ContourPlot('BestFit','ON'); close;
    sin2T4_contour  =  S.sin2T4_contour;
    mNu4Sq_contour = S.mNu4Sq_contour;
    sin2T4_contour_bf  =  S.sin2T4_bf;
    mNu4Sq_contour_bf = S.mNu4Sq_bf;
    %%
    
    chi2min = zeros(numel(mNuSq),1);
    chi2all = cell(numel(mNuSq),1);
    chi2all1 = cell(numel(mNuSq),1);
    mNu4Sq_iso = cell(numel(mNuSq),1);
    sin2T4_iso = cell(numel(mNuSq),1);
    mNu4Sq_bf  = zeros(numel(mNuSq),1);
    sin2T4_bf  = zeros(numel(mNuSq),1);
    
    for i=1:numel(mNuSq) 
        [M,c]= contour3(S.sin2T4,S.mNu4Sq,S.mNuSq,[mNuSq(i), mNuSq(i)],'LineWidth',2);  close;
        
        ExclIndex = find(M(1,:)==mNuSq(i));
        InclIndex = M(1,:)~=mNuSq(i);
        [X1,Y1] = meshgrid(S.mNu4Sq(1,:),S.sin2T4(:,1));
        % chi2min_tmp = zeros(numel(ExclIndex),1);
        % nur 1. contour nehmen
        
        %         for j=1:numel(ExclIndex)
        %             if j==numel(ExclIndex)
        %                 StopIdx = size(M,2);
        %             else
        %                 StopIdx = ExclIndex(j+1)-1;
        %             end
        %             chi2tmp = interp2(X1,Y1,S.chi2,M(2,ExclIndex(j):StopIdx),M(1,ExclIndex(j):StopIdx),'spline');
        %             chi2min_tmp(j) = min(chi2tmp);
        %
        %         end
        %        chi2min(i) = min(chi2min_tmp);
        mNu4Sq_iso{i}  = M(2,InclIndex);
        sin2T4_iso{i}  = M(1,InclIndex);
        
        chi2all{i} = interp2(X1,Y1,S.chi2,mNu4Sq_iso{i},sin2T4_iso{i},'spline');
        chi2all1 = interp2(X1,Y1,S.chi2,M(2,2:end),M(1,2:end),'spline');
        
        chi2min(i) =min(chi2all{i});
        BfIdx = find(chi2all{i}==min(chi2all{i}));
        
        mNutmp = mNu4Sq_iso{i};
        mNu4Sq_bf(i) = mNutmp(BfIdx);
        sin2T4tmp = sin2T4_iso{i};
        sin2T4_bf(i) = sin2T4tmp(BfIdx);
    end
    
    MakeDir(savedir);
    save(savefile,'mNuSq','chi2min','mNu4Sq_iso','sin2T4_iso','mNu4Sq_bf','sin2T4_bf',...
        'sin2T4_contour_bf','sin2T4_contour','mNu4Sq_contour_bf','mNu4Sq_contour');
    
end
%%
MinIdx = find(chi2min==min(chi2min));
chi2min_min = min(chi2min);
mNuSqErrUp  = interp1(chi2min(MinIdx:end),mNuSq(MinIdx:end),chi2min_min+1,'spline')-mNuSq(MinIdx);
mNuSqErrLow = mNuSq(MinIdx)-interp1(chi2min(1:MinIdx),mNuSq(1:MinIdx),chi2min_min+1,'spline');
mNuSqErr = mean([mNuSqErrUp,mNuSqErrLow]);
GetFigure;
if strcmp(DataType,'Real')
    mNuSqMin = mNuSq(MinIdx);
elseif strcmp(DataType,'Twin')
     mNuSqMin = 0;
end

[l,a]= boundedline(mNuSqMin.*ones(10,1),linspace(0,50,10),...
    [ones(10,1).*mNuSqErrLow,ones(10,1).*mNuSqErrUp],'orientation', 'horiz');
l.LineWidth = 2; l.Color = rgb('SteelBlue');
a.FaceColor = rgb('PowderBlue');
hold on;
%plot(mNuSq,min(chi2min)+1.*ones(numel(mNuSq),1),'--','Color',rgb('Black'),'LineWidth',1);
%plot(mNuSq,min(chi2min).*ones(numel(mNuSq),1),'--','Color',rgb('Black'),'LineWidth',1);
pchi2 = plot(mNuSq,chi2min,'-','LineWidth',3,'Color',rgb('Orange'),'MarkerSize',12);

xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'))
ylabel(sprintf('\\chi^2'))
PrettyFigureFormat('FontSize',22);

if strcmp(DataType,'Real')
xlim([-1 2.5])
ylim([min(chi2min)-1 max(chi2min)+1]);

leg = legend([pchi2,l],...
    sprintf('\\chi^2- profile data'),...
    sprintf('{\\itm}_\\nu^2 =  %.2f \\pm %.2f (^{-%.2f}_{+%.2f}) eV^2',mNuSq(MinIdx),mNuSqErr,mNuSqErrLow,mNuSqErrUp),...
    'Location','northwest');
elseif strcmp(DataType,'Twin')
  xlim([-1 1.5]) 
  ylim([0 5]);
  leg = legend([pchi2,l],...
    sprintf('\\chi^2- profile MC twin'),...
    sprintf('Sensitivity on {\\itm}_\\nu^2 = %.2f eV^2',mNuSqErrUp),...
    'Location','northwest');
end

PrettyLegendFormat(leg);


if strcmp(SavePlt,'ON')
    pltdir =  [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/plots/'];
    MakeDir(pltdir);
    pltname = strrep(strrep(savefile,'results','plots'),'.mat','.png');
    print(pltname,'-dpng','-r350');
    fprintf('save plot to %s \n',pltname)
end






% combine results from randomized grid search: ksn1+ksn2

% create smaller file for result of grids searches on randomized data

InterpMode = 'lin';
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
RecomputeFlag = 'OFF';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn21_WilksTheorem_NullHypothesis_Interp%s.mat',...
        savedir,InterpMode);
else
    savefile = sprintf('%sksn21_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode);
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    fprintf('savefile already created \n');
   d = importdata(savefile);
   load(savefile)
else
    
    %% configure RunAnalysis object  + SterileAnalysis object
    DataType = 'Twin';
    freePar = 'E0 Norm Bkg';
    nGridSteps = 25;
    range = 40;
    chi2 = 'chi2CMShape';
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
    % configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'range',range,...
        'RandMC','OFF',...
        'Twin_mNu4Sq',Twin_mNu4Sq,...
        'Twin_sin2T4',Twin_sin2T4};
    
    S = SterileAnalysis(SterileArg{:});
    
    %% KSN-2
    %% init variables: ksn2 first part
    randMC =[1001:1260,1294:1300,1349:1500];
    NrandMC1 = numel(randMC);
    chi2_grid      = cell(NrandMC1,1);
    TBDIS_mc      = zeros(NrandMC1,A.ModelObj.nqU);
    chi2_nullk2   = zeros(NrandMC1,1);
    
    S.InterpMode = InterpMode;
    S.LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
    
    % first part: load files
    progressbar('Merge files for WT');
    for i=1:NrandMC1
        progressbar(i/NrandMC1);
        S.RandMC= randMC(i);
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Minm4Sq',1,'Maxm4Sq',38^2)
        chi2_grid{i}   = S.chi2;
        TBDIS_mc(i,:) = S.RandMC_TBDIS;
        S.RandMC_TBDIS = [];
        chi2_nullk2(i) = S.chi2_Null;
    end
    mNu4Sq_grid = S.mNu4Sq;
    sin2T4_grid = S.sin2T4;
    dof = S.dof;
    
    % remove duplicants from first part
    % find grid with idential TBDIS_mc
    [~,IdxNonId] = unique(TBDIS_mc(1:numel(randMC),:),'rows');
    IdxNonId = sort(IdxNonId);
    chi2_grid  = chi2_grid(IdxNonId);
    TBDIS_mc   = TBDIS_mc(IdxNonId,:);
    chi2_nullk2 = chi2_nullk2(IdxNonId);
    
    %% ksn-2 second part
    randMC_new  = 1:1250;
    NrandMC2 = numel(randMC_new);
    TBDIS_mc2      = zeros(NrandMC2,A.ModelObj.nqU);
    chi2_grid2     = cell(NrandMC2,1);
    chi2_nullk22     = zeros(NrandMC2,1);
    
    dspectra = ksn2_GetRandTritiumSpectrum('H0');
    progressbar('Merge files for WT - part 2');
    for i=1:NrandMC2
        progressbar(i/NrandMC2);
        S.RandMC= randMC_new(i);
        S.RandMC_TBDIS = dspectra.TBDIS_mc(:,randMC_new(i));
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Minm4Sq',1,'Maxm4Sq',38^2)
        chi2_grid2{i}  = S.chi2;
        TBDIS_mc2(i,:) = S.RandMC_TBDIS;
        S.RandMC_TBDIS = [];
        chi2_nullk22(i) = S.chi2_Null;
    end
    mNu4Sq_grid2 = S.mNu4Sq;
    sin2T4_grid2 = S.sin2T4;
    %% KSN2 merge 1st and 2nd part
    chi2_grid         = [chi2_grid;chi2_grid2];
    TBDIS_mc          = [TBDIS_mc;TBDIS_mc2];
    chi2_nullk2         = [chi2_nullk2;chi2_nullk22];
    %% check if they have same binning
    if sum(sum(sin2T4_grid~=sin2T4_grid2))>0 || sum(sum(mNu4Sq_grid~=mNu4Sq_grid2))>0
        fprintf(2, 'KSN-1 and KSN-2 do not have  the same binning - return \n');
        return;
    end
    %% load ksn1
    RandMC = 1:1500;
    chi2_gridksn1     = cell(numel(RandMC),1);
    chi2_nullksn1 = zeros(numel(RandMC),1);
    for i=RandMC
        progressbar(i/numel(RandMC))
        [mnu4Sq_tmp,sin2T4_tmp,chi2_tmp,~,~,FitResults_Null] = KSN1GridSearch('range',range,...
            'nGridSteps',nGridSteps,...
            'chi2','chi2CMShape',...
            'DataType','Twin',...
            'freePar','E0 Bkg Norm',...
            'RunList','KNM1',...
            'SmartGrid','OFF',...
            'RecomputeFlag','OFF',...
            'RandMC',i,...
            'SysBudget',24,...
            'Twin_mNu4Sq',Twin_mNu4Sq,...
            'Twin_sin2T4',Twin_sin2T4);
        S.chi2 = chi2_tmp;
        S.mNu4Sq = mnu4Sq_tmp;
        S.sin2T4 = sin2T4_tmp;
        S.mNuSq = zeros(nGridSteps,nGridSteps);
        S.E0 = zeros(nGridSteps,nGridSteps);
        S.Interp1Grid('Minm4Sq',1,'Maxm4Sq',38^2)
        chi2_gridksn1{i}  = S.chi2;
        chi2_nullksn1(i) = FitResults_Null.chi2min;
    end
    sin2T4_gridksn1 = S.sin2T4;
    mNu4Sq_gridksn1 = S.mNu4Sq;
    
    if sum(sum(sin2T4_grid~=sin2T4_gridksn1))>0 || sum(sum(mNu4Sq_grid~=mNu4Sq_gridksn1))>0
        fprintf(2, 'KSN-1 and KSN-2 do not have  the same binning - return \n');
        return;
    end
    %% combine:
    nCombi = 1500;
    
    mNu4Sq_bf      = zeros(nCombi,1);
    sin2T4_bf      = zeros(nCombi,1);
    chi2_bf        = zeros(nCombi,1);
    chi2_null      = zeros(nCombi,1);
    mNu4Sq_contour = cell(nCombi,1);
    sin2T4_contour = cell(nCombi,1);
    
    for i=1:nCombi
        chi2sum = chi2_grid{i}+chi2_gridksn1{i};
        chi2_bf(i) = min(min(chi2sum));
        DeltaChi2sum = chi2sum-chi2_bf(i);
        S.sin2T4 = sin2T4_grid;
        S.mNu4Sq = mNu4Sq_grid;
        S.chi2 = DeltaChi2sum;
        S.chi2_ref = 0;
        chi2_null(i)  = chi2_nullk2(i) + chi2_nullksn1(i);
        
        S.ContourPlot('BestFit','ON');close;
        mNu4Sq_bf(i) = S.mNu4Sq_bf;
        sin2T4_bf(i) = S.sin2T4_bf;
        mNu4Sq_contour{i} = S.mNu4Sq_contour;
        sin2T4_contour{i} = S.sin2T4_contour;
    end
    
    chi2_delta = chi2_null-chi2_bf;
     ClosedLog95 =  chi2_delta>= GetDeltaChi2(0.95,2);
    save(savefile,'mNu4Sq_contour','sin2T4_contour',...
        'chi2_bf','chi2_null','chi2_delta',...
        'mNu4Sq_bf','sin2T4_bf','ClosedLog95');
end
%% interpolate



GetFigure
for i=1:100
    
    sin2T4_tmp = cell2mat(sin2T4_contour(i));
    mNu4Sq_tmp = cell2mat(mNu4Sq_contour(i));
%     if size(sin2T4_tmp
    
    plot(sin2T4_tmp,mNu4Sq_tmp,':');
    hold on;
end
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off
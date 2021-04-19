% create smaller file for result of grids searches on randomized data
RecomputeFlag = 'ON';
randMC = 1:1e3;
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,numel(randMC));
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,numel(randMC));
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    fprintf('savefile already created \n');
    d = importdata(savefile);
else
    %% configure RunAnalysis object
    chi2 = 'chi2CMShape';
    DataType = 'Twin';
    freePar = 'E0 Norm Bkg';
    nGridSteps = 25;
    range = 40;
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
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'range',range,...
        'RandMC','OFF',...
        'Twin_mNu4Sq',Twin_mNu4Sq,...
        'Twin_sin2T4',Twin_sin2T4};
    
    S = SterileAnalysis(SterileArg{:});
    %% init variables
    mNu4Sq_bf      = zeros(numel(randMC),1);
    sin2T4_bf      = zeros(numel(randMC),1);
    chi2_bf        = zeros(numel(randMC),1);
    chi2_null      = zeros(numel(randMC),1);
    chi2_delta     = zeros(numel(randMC),1);
    mNu4Sq_contour = cell(numel(randMC),1);
    sin2T4_contour = cell(numel(randMC),1);
    
    S.InterpMode = 'lin';
    S.LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
    
    % load files
    for i=1:numel(randMC)
        progressbar(i/numel(randMC));
        S.RandMC= randMC(i);
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('Maxm4Sq',38^2);
        S.ContourPlot; close;
       
        S.FindBestFit('Mode','Def');
        S.FindBestFit('Mode','Imp');
        
        mNu4Sq_bf(i) = S.mNu4Sq_bf;
        sin2T4_bf(i) = S.sin2T4_bf;
        chi2_bf(i)   = S.chi2_bf;
        chi2_null(i) = S.chi2_Null;
        chi2_delta(i) = S.chi2_Null-S.chi2_bf;
        
        mNu4Sq_contour{i} = S.mNu4Sq_contour;
        sin2T4_contour{i} = S.sin2T4_contour;
    end
    
    dof = S.dof;
    
    %% also load expected contour from Asimov Twin
    S.InterpMode = 'spline';
    S.RandMC= 'OFF';
    S.nGridSteps = 30;
    S.LoadGridFile('mNu4SqTestGrid',5,'ExtmNu4Sq','ON');
    S.Interp1Grid;
    S.ContourPlot; close;
    mNu4Sq_contour_Asimov = S.mNu4Sq_contour;
    sin2T4_contour_Asimov = S.sin2T4_contour;
    
    MakeDir(savedir);
    save(savefile,'mNu4Sq_bf','sin2T4_bf','chi2_bf','chi2_null',...
        'chi2_delta','mNu4Sq_contour','sin2T4_contour',...
        'mNu4Sq_contour_Asimov','sin2T4_contour_Asimov');
    fprintf('save file to %s \n',savefile);
    d = importdata(savefile);
end

%% to some further calculations
if ~isfield(d,'ClosedLog95')
    % number of significant best fits
    nContours = numel(randMC);
    ClosedLog95 =  chi2_delta>= GetDeltaChi2(0.95,2);
    ClosedFrac95 = sum(ClosedLog95)/nContours;
    RelErr = @(n,p) sqrt((n*p*(1-p)))/n;
    ClosedFrac95RelErr = RelErr(nContours,ClosedFrac95);
    fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  +- %.1f  (%.0f out of %.0f)\n',...
        95,ClosedFrac95*100,ClosedFrac95RelErr*100,sum(ClosedLog95),nContours);
    save(savefile,'ClosedLog95','ClosedFrac95','ClosedFrac95RelErr','-append');
end


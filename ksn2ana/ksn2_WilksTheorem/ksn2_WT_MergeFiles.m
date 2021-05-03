% create smaller file for result of grids searches on randomized data
RecomputeFlag = 'ON';
Hypothesis = 'H0';
InterpMode = 'lin';

switch Hypothesis
    case 'H0'
        %randMC = [1:1e3];
       % randMC =[1001:1258,1296:1300,1351:1500];%11:1e3;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1' 
        randMC = 1:1e3;
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples.mat',savedir,InterpMode,numel(randMC));
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC));
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    fprintf('savefile already created \n');
    d = importdata(savefile);
else
    %% configure RunAnalysis object 
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
    ReCalc_chi2Null   = zeros(numel(randMC),1); 
    ReCalc_chi2Null_i = zeros(numel(randMC),1);  
    chi2_delta     = zeros(numel(randMC),1);
    mNu4Sq_contour = cell(numel(randMC),1);
    sin2T4_contour = cell(numel(randMC),1);
    TBDIS_mc      = zeros(numel(randMC),A.ModelObj.nqU);
     
    S.InterpMode = InterpMode;
    S.LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
    
    % load files
    progressbar('Merge files for WT');
    for i=1:numel(randMC)
        progressbar(i/numel(randMC));
        S.RandMC= randMC(i);
        S.LoadGridFile(S.LoadGridArg{:});
        
        if strcmp(InterpMode,'Mix')
            S.Interp1Grid('Maxm4Sq',36^2);% In "Mix" Mode -> Maxm4Sq is threshold for spline/linear interpolation
        else
            S.Interp1Grid;
        end
        
        if strcmp(Hypothesis,'H1')
            S.FindBestFit;
            if S.mNu4Sq_bf>1000
                InterpMode_i = S.InterpMode;
                S.InterpMode = 'spline';
                S.LoadGridFile(S.LoadGridArg{:});
                S.Interp1Grid('Maxm4Sq',500,'Minm4Sq',1);
                S.InterpMode = InterpMode_i;
            end
        end
        
        S.ContourPlot('BestFit','ON'); close;
        mNu4Sq_contour{i} = S.mNu4Sq_contour;
        sin2T4_contour{i} = S.sin2T4_contour;
        
        %    S.FindBestFit('Mode','Imp');
        mNu4Sq_bf(i) = S.mNu4Sq_bf;
        sin2T4_bf(i) = S.sin2T4_bf;
        chi2_bf(i)   = S.chi2_bf;
        chi2_null(i) = S.chi2_Null;
        chi2_delta(i) = chi2_null(i)-S.chi2_bf;
        TBDIS_mc(i,:) = S.RandMC_TBDIS;
        
        S.RandMC_TBDIS = [];
        
       if  randMC(i)<=1e3
        df = S.GridFilename(S.LoadGridArg{:});
        d = importdata(df);
        ReCalc_chi2Null_i(i) = d.ReCalc_chi2Null_i;
        ReCalc_chi2Null(i) = d.ReCalc_chi2Null;
       end
    end
    
    dof = S.dof;
    MakeDir(savedir);
    
        
    % also load expected contour from Asimov Twin
    S.InterpMode = 'Mix';%'spline';
    S.RandMC= 'OFF';
    S.nGridSteps = 30;
    
    S.LoadGridFile('mNu4SqTestGrid',5,'ExtmNu4Sq','ON');
    S.Interp1Grid;
    S.ContourPlot; close;
    
    mNu4Sq_contour_Asimov = S.mNu4Sq_contour;
    sin2T4_contour_Asimov = S.sin2T4_contour;
    if max(randMC)<=1e3
    save(savefile,'mNu4Sq_bf','sin2T4_bf','chi2_bf','chi2_null',...
        'ReCalc_chi2Null','ReCalc_chi2Null_i',...
        'chi2_delta','mNu4Sq_contour','sin2T4_contour',....
        'mNu4Sq_contour_Asimov','sin2T4_contour_Asimov','TBDIS_mc');
    else
          
        save(savefile,'mNu4Sq_bf','sin2T4_bf','chi2_bf','chi2_null',...
        'chi2_delta','mNu4Sq_contour','sin2T4_contour',....
        'mNu4Sq_contour_Asimov','sin2T4_contour_Asimov','TBDIS_mc');
    end
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

%% find grid with idential TBDIS_mc
if ~isfield(d,'IdxNonId')
    [TBDIS_mcNonId,IdxNonId] = unique(TBDIS_mc,'rows');
    ClosedLog95_NonId =  chi2_delta(IdxNonId)>= GetDeltaChi2(0.95,2);
    ClosedFrac95_NonId = sum(ClosedLog95_NonId)/numel(IdxNonId);
    RelErr = @(n,p) sqrt((n*p*(1-p)))/n;
    ClosedFrac95RelErr_NonId = RelErr(numel(IdxNonId),ClosedFrac95_NonId);
    fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  +- %.1f  (%.0f out of %.0f)\n',...
        95,ClosedFrac95_NonId*100,ClosedFrac95RelErr_NonId*100,sum(ClosedLog95_NonId),numel(IdxNonId));
    save(savefile,'IdxNonId','-append');
end
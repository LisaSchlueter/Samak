% create smaller file for result of grids searches on randomized data
RecomputeFlag = 'OFF';
Hypothesis = 'H2';
InterpMode = 'Mix';%'lin';
chi2 = 'chi2CMShape';
randMC = [1:1500];
NrandMC = numel(randMC);
if strcmp(Hypothesis,'H1')
    Twin_sin2T4 = 0.0240;
    Twin_mNu4Sq = 92.7;
    mNu4SqTestGrid = 2;
else
    Twin_sin2T4 = 0.07;
    Twin_mNu4Sq = 20;
    mNu4SqTestGrid = 'OFF';
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
    savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,NrandMC);

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
    
    %% init variables: first part
    mNu4Sq_bf      = zeros(NrandMC,1);
    sin2T4_bf      = zeros(NrandMC,1);
    chi2_bf        = zeros(NrandMC,1);
    chi2_null      = zeros(NrandMC,1);  
    chi2_delta     = zeros(NrandMC,1);
    mNu4Sq_contour = cell(NrandMC,1);
    sin2T4_contour = cell(NrandMC,1);
    TBDIS_mc      = zeros(NrandMC,A.ModelObj.nqU);
    
    S.InterpMode = InterpMode;
 
    
    dspectra = ksn2_GetRandTritiumSpectrum(Hypothesis);
    progressbar('Merge files for WT');
    for i=1:NrandMC
        progressbar(i/NrandMC);
        S.RandMC= randMC(i);
        if randMC(i)>1000
            S.LoadGridArg = {'mNu4SqTestGrid','OFF','ExtmNu4Sq','ON'};
        else
            S.LoadGridArg = {'mNu4SqTestGrid',mNu4SqTestGrid,'ExtmNu4Sq','ON'};
        end
        S.RandMC_TBDIS = dspectra.TBDIS_mc(:,randMC(i));
        S.LoadGridFile(S.LoadGridArg{:});
        if strcmp(InterpMode,'Mix')
            S.Interp1Grid('Maxm4Sq',700);% In "Mix" Mode -> Maxm4Sq is threshold for spline/linear interpolation
        else
            S.Interp1Grid('Maxm4Sq',36^2);
        end
        
        S.ContourPlot('BestFit','ON'); close;
        mNu4Sq_contour{i} = S.mNu4Sq_contour;
        sin2T4_contour{i} = S.sin2T4_contour;
        
        mNu4Sq_bf(i) = S.mNu4Sq_bf;
        sin2T4_bf(i) = S.sin2T4_bf;
        chi2_bf(i)   = S.chi2_bf;
        chi2_null(i) = S.chi2_Null;
        chi2_delta(i) = chi2_null(i)-S.chi2_bf;
        TBDIS_mc(i,:) = S.RandMC_TBDIS;
        
        S.RandMC_TBDIS = [];
    end
    % also load expected contour from Asimov Twin
   
    S.RandMC= 'OFF';
    S.nGridSteps = 30;
    if strcmp(Hypothesis,'H1')
         S.InterpMode = 'Mix';%'spline';
        S.LoadGridFile('mNu4SqTestGrid',5,'ExtmNu4Sq','ON');
         S.Interp1Grid;
    elseif  strcmp(Hypothesis,'H2')
         S.InterpMode = 'Mix';
        S.LoadGridFile('mNu4SqTestGrid','OFF','ExtmNu4Sq','ON');
         S.Interp1Grid('Maxm4Sq',30^2);
    end
   
    S.ContourPlot; close;
    
    mNu4Sq_contour_Asimov = S.mNu4Sq_contour;
    sin2T4_contour_Asimov =  S.sin2T4_contour;
    
    if ismember(Hypothesis,{'H1','H2'})
        S.FindBestFit;
        mNu4Sq_contour_Asimov_bf = S.mNu4Sq_bf;
        sin2T4_contour_Asimov_bf =  S.sin2T4_bf;
    end
    % save
     MakeDir(savedir);
    save(savefile,'mNu4Sq_bf','sin2T4_bf','chi2_bf','chi2_null',...
        'chi2_delta','mNu4Sq_contour','sin2T4_contour',....
        'mNu4Sq_contour_Asimov','sin2T4_contour_Asimov','TBDIS_mc');
    
    if ismember(Hypothesis,{'H1','H2'}) 
        save(savefile,'mNu4Sq_contour_Asimov_bf','sin2T4_contour_Asimov_bf','-append');
    end
    
    fprintf('save file to %s \n',savefile);
    d = importdata(savefile);
    
end

%% to some further calculations
if ~isfield(d,'ClosedLog95')
    % number of significant best fits
    nContours = numel(mNu4Sq_bf);
    ClosedLog95 =  chi2_delta>= GetDeltaChi2(0.95,2);
    ClosedFrac95 = sum(ClosedLog95)/nContours;
    RelErr = @(n,p) sqrt((n*p*(1-p)))/n;
    ClosedFrac95RelErr = RelErr(nContours,ClosedFrac95);
    fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  +- %.1f  (%.0f out of %.0f)\n',...
        95,ClosedFrac95*100,ClosedFrac95RelErr*100,sum(ClosedLog95),nContours);
    save(savefile,'ClosedLog95','ClosedFrac95','ClosedFrac95RelErr','-append');
end

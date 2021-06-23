% create smaller file for result of grids searches on randomized data
RecomputeFlag = 'OFF';
InterpMode = 'spline';
MergeNew = 'ON';
RmDuplicates = 'ON';
randMC =[1001:1260,1294:1300,1349:1500];
randMC_new  = 1:1250;
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
chi2 = 'chi2CMShape';

%%
if strcmp(MergeNew,'ON')
    MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
      NrandMC = numel(randMC)+numel(randMC_new);
else
    MergeStr = '';
    NrandMC = numel(randMC);
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,numel(randMC),MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,numel(randMC),MergeStr,RmDuplicates);
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
    
    %% init variables: first part
    NrandMC1 = numel(randMC);
    mNu4Sq_bf      = zeros(NrandMC1,1);
    sin2T4_bf      = zeros(NrandMC1,1);
    chi2_bf        = zeros(NrandMC1,1);
    chi2_null      = zeros(NrandMC1,1);  
    chi2_delta     = zeros(NrandMC1,1);
    mNu4Sq_contour = cell(NrandMC1,1);
    sin2T4_contour = cell(NrandMC1,1);
    TBDIS_mc      = zeros(NrandMC1,A.ModelObj.nqU);
    
    S.InterpMode = InterpMode;
    S.LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
    
    savefileOld = [extractBefore(savefile,'_MergeNew'),'.mat'];
    if exist(savefileOld,'file')
        load(savefileOld)
        fprintf('load %s \n',savefileOld);
    else
        
        % first part: load files
        progressbar('Merge files for WT');
        for i=1:NrandMC1
            progressbar(i/NrandMC1);
            S.RandMC= randMC(i);
            S.LoadGridFile(S.LoadGridArg{:});
            
            if strcmp(InterpMode,'Mix')
                S.Interp1Grid('Maxm4Sq',36^2);% In "Mix" Mode -> Maxm4Sq is threshold for spline/linear interpolation
            else
                S.Interp1Grid;
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
        
        dof = S.dof;
        save(savefileOld,'mNu4Sq_contour','sin2T4_contour','mNu4Sq_bf','chi2_bf','chi2_null','chi2_delta','TBDIS_mc')
        
    end
    %% remove duplicants from first part
    % find grid with idential TBDIS_mc
    if strcmp(RmDuplicates,'ON')
        [~,IdxNonId] = unique(TBDIS_mc(1:numel(randMC),:),'rows');
        IdxNonId = sort(IdxNonId);
        mNu4Sq_contour = mNu4Sq_contour(IdxNonId);
        sin2T4_contour = sin2T4_contour(IdxNonId);
        
        mNu4Sq_bf     = mNu4Sq_bf(IdxNonId);
        sin2T4_bf     = sin2T4_bf(IdxNonId);
        chi2_bf       = chi2_bf(IdxNonId);
        chi2_null     = chi2_null(IdxNonId);
        chi2_delta    = chi2_delta(IdxNonId);
        TBDIS_mc      = TBDIS_mc(IdxNonId,:);  
    end
    
    
    %% second part
    if strcmp(MergeNew,'ON')
        
        NrandMC2 = numel(randMC_new);
        mNu4Sq_bf2      = zeros(NrandMC2,1);
        sin2T4_bf2      = zeros(NrandMC2,1);
        chi2_bf2        = zeros(NrandMC2,1);
        chi2_null2      = zeros(NrandMC2,1);
        chi2_delta2     = zeros(NrandMC2,1);
        mNu4Sq_contour2 = cell(NrandMC2,1);
        sin2T4_contour2 = cell(NrandMC2,1);
        TBDIS_mc2      = zeros(NrandMC2,A.ModelObj.nqU);

        dspectra = ksn2_GetRandTritiumSpectrum('H0');
        progressbar('Merge files for WT - part 2');   
        for i=1:NrandMC2
            progressbar(i/NrandMC2);
            S.RandMC= randMC_new(i);
            S.RandMC_TBDIS = dspectra.TBDIS_mc(:,randMC_new(i));
            S.LoadGridFile(S.LoadGridArg{:});
            if strcmp(InterpMode,'Mix')
                S.Interp1Grid('Maxm4Sq',36^2);% In "Mix" Mode -> Maxm4Sq is threshold for spline/linear interpolation
            else
                S.Interp1Grid;
            end
            
            S.ContourPlot('BestFit','ON'); close;
            mNu4Sq_contour2{i} = S.mNu4Sq_contour;
            sin2T4_contour2{i} = S.sin2T4_contour;
            
            mNu4Sq_bf2(i) = S.mNu4Sq_bf;
            sin2T4_bf2(i) = S.sin2T4_bf;
            chi2_bf2(i)   = S.chi2_bf;
            chi2_null2(i) = S.chi2_Null;
            chi2_delta2(i) = chi2_null2(i)-S.chi2_bf;
            TBDIS_mc2(i,:) = S.RandMC_TBDIS;
            
            S.RandMC_TBDIS = [];
        end  
            %% merge 1st and 2nd part
     
            mNu4Sq_bf       = [mNu4Sq_bf;mNu4Sq_bf2];
            sin2T4_bf       = [sin2T4_bf;sin2T4_bf2];
            chi2_bf         = [chi2_bf;chi2_bf2];
            chi2_null       = [chi2_null;chi2_null2];
            chi2_delta      = [chi2_delta;chi2_delta2];
            TBDIS_mc        = [TBDIS_mc;TBDIS_mc2];
            mNu4Sq_contour  = [mNu4Sq_contour;mNu4Sq_contour2];
            sin2T4_contour = [sin2T4_contour;sin2T4_contour2];
            
          
    end   
    %% also load expected contour from Asimov Twin
    if strcmp(InterpMode,'spline')
        S.InterpMode = 'spline';
    else
        S.InterpMode = 'Mix';%'spline';
    end
    S.RandMC= 'OFF';
    S.nGridSteps = 30;
    
    S.LoadGridFile('mNu4SqTestGrid',5,'ExtmNu4Sq','ON');
    S.Interp1Grid;
    S.ContourPlot; close;
    
    mNu4Sq_contour_Asimov = S.mNu4Sq_contour;
    sin2T4_contour_Asimov = S.sin2T4_contour;
    
    %% save
     MakeDir(savedir);
    save(savefile,'mNu4Sq_bf','sin2T4_bf','chi2_bf','chi2_null',...
        'chi2_delta','mNu4Sq_contour','sin2T4_contour',....
        'mNu4Sq_contour_Asimov','sin2T4_contour_Asimov','TBDIS_mc');
    
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

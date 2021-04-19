% do combi fit knm1+2 for sterile neutrinos
% with new ksn1 model 
% does only make sense for data (other twins have to be re-calculated)
RecomputeFlag = 'OFF';

%% KNM-1 Model
chi2          = 'chi2CMShape';
DataType      = 'Real';
nGridSteps    = 10;
freePar       = 'mNu E0 Bkg Norm';
range         = 40;
savedir = [getenv('SamakPath'),'/SterileAnalysis/GridSearchFiles/Combi/',DataType,'/'];
MakeDir(savedir)
savename = sprintf('%sKSN12Combi_ReAna_GridSearch_%s_%s_Uniform_%s_%.0fnGrid.mat',savedir,DataType,strrep(freePar,' ',''),chi2,nGridSteps);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
    tStart = tic;
    FSDFlag           = 'KNM2_0p5eV';
    RunList           = 'KNM1';
    DopplerEffectFlag = 'FSD';
    PullFlag          = 99;
    
    K1 = MultiRunAnalysis('RunList',RunList,...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',200,...
        'AngularTFFlag','ON',...
        'DopplerEffectFlag',DopplerEffectFlag,...
        'BKG_PtSlope',-2.2.*1e-06,...
        'PullFlag',PullFlag);
    
    K1.exclDataStart = K1.GetexclDataStart(range);
    K1.i_qUOffset    = zeros(1,K1.ModelObj.nPixels);
    K1.i_mTSq        = zeros(1,K1.ModelObj.nPixels);
    
    %% KNM-2 Model
    AnaFlag               = 'StackPixel';
    BKG_PtSlope           = 3*1e-06;
    TwinBias_BKG_PtSlope  = 3*1e-06;
    SysBudget             = 40;
    NonPoissonScaleFactor = 1.112;
    SigmaSq               =  0.0124+0.0025;
    
    % Init Model Object and covariance matrix object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',PullFlag,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    K2 = MultiRunAnalysis(RunAnaArg{:});
    K2.exclDataStart = K2.GetexclDataStart(range);
    K2.i_qUOffset = zeros(1,K2.ModelObj.nPixels);
    K2.i_mTSq = zeros(1,K2.ModelObj.nPixels);
    
    %% Combine
    Data = [K1.RunData.qU(2:end),K2.RunData.qU,...
        K1.RunData.TBDIS(2:end),K2.RunData.TBDIS,...
        K1.RunData.TBDISE(2:end),K2.RunData.TBDISE];
    
%% cov mats
    K1.ModelObj.nPixels = 2; % to make FITC read data properly
    if strcmp(K1.chi2,'chi2Stat')
        if ~contains(K1.fixPar,'fix 3 ;')
            K1.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
        end
        
        [StatCM, StatCMFrac] = K1.ComputeCM_StatPNP(varargin);
        K1.FitCM = StatCM;
        K1.FitCMShape = StatCM;
        K1.FitCMFrac = StatCMFrac;
        K1.FitCMFracShape = StatCMFrac;
        
        if ~contains(K2.fixPar,'fix 3 ;')
            K2.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
        end
        
        [StatCM2, StatCMFrac2] = K2.ComputeCM_StatPNP(varargin);
        K2.FitCM      = StatCM2;
        K2.FitCMShape = StatCM2;
        K2.FitCMFrac = StatCMFrac2;
        K2.FitCMFracShape = StatCMFrac2;
        
    end
    
   % add knm1+2 cov mats 
    nqU = size(Data,1);
    COVMAT          = diag(zeros(2*nqU,2*nqU));
    COVMATShape     = diag(zeros(2*nqU,2*nqU));
    COVMATFracShape = diag(zeros(2*nqU,2*nqU));
    
    COVMAT(1:nqU,1:nqU) = K1.FitCM(2:end,2:end);
    COVMAT(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCM;
    COVMATShape(1:nqU,1:nqU) = K1.FitCMShape(2:end,2:end);
    COVMATShape(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCMShape;
    COVMATFracShape(1:nqU,1:nqU) = K1.FitCMFracShape(2:end,2:end);
    COVMATFracShape(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCMFracShape;
    
%% prepare fit
    FitCArg = {'SO',K1.ModelObj,...
        'SO2',K2.ModelObj,...
        'DATA',Data,'fitter',K1.fitter,...
        'chi2name',K1.chi2,'minuitOpt',K1.minuitOpt,...
        'COVMAT', real(COVMAT),'COVMATFrac', real(K1.FitCMFrac),...
        'COVMATShape', real(COVMATShape),'COVMATNorm',K1.FitCMNorm,...
        'COVMATFracShape',real(COVMATFracShape),...
        'pulls',K1.pulls,...
        'pullFlag',K1.pullFlag,...
        'fixPar','',...
        'exclDataStart',K1.exclDataStart,...
        'exclDataStop',K1.exclDataStop,...
        'i_mnu',K1.i_mnu,...
        'i_Q',K1.i_Q,...
        'i_Q2',K2.i_Q,...
        'i_B',K1.i_B,...
        'i_N',K1.i_N,...
        'i_DTGS',K1.i_DTGS,...
        'i_DTES',K1.i_DTES,...
        'i_HTGS',K1.i_HTGS,...
        'i_HTES',K1.i_HTES,...
        'i_TTGS',K1.i_TTGS,...
        'i_TTES',K1.i_TTES,...
        'i_qUOffset',K1.i_qUOffset,...
        'i_mTSq',K1.i_mTSq,...
        'i_FracTm',0};
    
    %% reference fit with m4=0 sint4Sq =0 (Null hypothesis)
   F = FITC(FitCArg{:});
  FitResult_Null = struct(...
      'par',F.RESULTS{1},....
      'err',F.RESULTS{2},....
      'chi2min',F.RESULTS{3},...
      'errmat',F.RESULTS{4},...
      'dof',F.RESULTS{5});
  
  % define grid
  sin2T4Max = 0.5;
  if nGridSteps<4 % for testing
      mnu4Sq = logspace(0,log10((range)^2),nGridSteps)';
  else
      mNu4Sq_ex = [0.1;0.35;0.7];
      nGridSteps = nGridSteps-3;
      mnu4SqSmall  = logspace(0,log10((range-11)^2),nGridSteps-10)';
      mnu4SqLarge  = logspace(log10((range-10)^2),log10((range)^2),10)';
      mnu4Sq       = sort([mNu4Sq_ex;mnu4SqSmall;mnu4SqLarge]);
  end

sin2T4      = logspace(-3,log10(sin2T4Max),nGridSteps);
mnu4Sq      = repmat(mnu4Sq,1,nGridSteps);
sin2T4      = repmat(sin2T4,nGridSteps,1);

%% make copy of model for parallel computing
D1 = copy(repmat(K1,nGridSteps^2,1));
D2 = copy(repmat(K2,nGridSteps^2,1));

D1 = reshape(D1,numel(D1),1);
D2 = reshape(D2,numel(D2),1);
%% scan over msq4-sin2t4 grid
chi2Grid       = zeros(nGridSteps^2,1);
FitResultsGrid = cell(nGridSteps^2,1);
mnu4Sq_Grid    = reshape(mnu4Sq',nGridSteps^2,1);
sin2T4_Grid    = reshape(sin2T4',nGridSteps^2,1);

    fixPartmp = K1.fixPar;
    DataTypetmp = DataType;
    
    parfor i= 1:(nGridSteps^2)
       
        D1(i).SimulateStackRuns;
        D2(i).SimulateStackRuns;
        D1(i).ModelObj.nPixels = 2; %
        
        % set sterile parameter
        D1(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
        D2(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
        
        % fit
        F = FITC('SO',D1(i).ModelObj,...
        'SO2',D2(i).ModelObj,...
        'DATA',Data,'fitter',D1(i).fitter,...
        'chi2name',D1(i).chi2,'minuitOpt',D1(i).minuitOpt,...
        'COVMAT', real(COVMAT),'COVMATFrac', real(D1(i).FitCMFrac),...
        'COVMATShape', real(COVMATShape),'COVMATNorm',D1(i).FitCMNorm,...
        'COVMATFracShape',real(COVMATFracShape),...
        'pulls',D1(i).pulls,...
        'pullFlag',D1(i).pullFlag,...
        'fixPar','',...
        'exclDataStart',D1(i).exclDataStart,...
        'exclDataStop',D1(i).exclDataStop,...
        'i_mnu',D1(i).i_mnu,...
        'i_Q',D1(i).i_Q,...
        'i_Q2',D2(i).i_Q,...
        'i_B',D1(i).i_B,...
        'i_N',D1(i).i_N,...
        'i_DTGS',D1(i).i_DTGS,...
        'i_DTES',D1(i).i_DTES,...
        'i_HTGS',D1(i).i_HTGS,...
        'i_HTES',D1(i).i_HTES,...
        'i_TTGS',D1(i).i_TTGS,...
        'i_TTES',D1(i).i_TTES,...
        'i_qUOffset',D1(i).i_qUOffset,...
        'i_mTSq',D1(i).i_mTSq,...
        'i_FracTm',0);
        
        FitResult = struct(...
            'par',F.RESULTS{1},....
            'err',F.RESULTS{2},....
            'chi2min',F.RESULTS{3},...
            'errmat',F.RESULTS{4},...
            'dof',F.RESULTS{5});
        
        chi2Grid(i) = F.RESULTS{3};
        FitResultsGrid{i} = FitResult;
    end
    
    chi2       = reshape(chi2Grid,nGridSteps,nGridSteps);
    FitResults = reshape(FitResultsGrid,nGridSteps,nGridSteps);
    
    if min(min(chi2))<FitResult_Null.chi2min
        chi2_ref = min(min(chi2));
    else
        chi2_ref = FitResult_Null.chi2min;
    end
    
    tCpuHour = (cputime-tStart)/60; % cpu time in hours
    %% save
    save(savename,'chi2_ref',...%'FitResults_ref'
        'chi2','mnu4Sq','sin2T4','FitResults','FitResults_Null','tCpuHour');
    fprintf('save file to %s \n',savename);

end


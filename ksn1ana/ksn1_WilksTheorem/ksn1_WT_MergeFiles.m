% contour with randomized twins
%% settings
range = 40;%
chi2Str = 'chi2CMShape';
InterpMode = 'lin';
savedir = [getenv('SamakParh'),'ksn1ana/ksn1_WilksTheorem/results/'];
MakeDir(savedir);
savename = sprintf('%sksn1_WilksTheorem_%.0frange_%s_%s.mat',savedir,range,chi2Str,InterpMode);

if exist(savename,'file')
    d = importdata(savename);
else
    
    CL = 95;
    
    nGridSteps = 25;
    
    DataType = 'Twin';
    freePar = 'E0 Bkg Norm';
    RunList = 'KNM1';
    SmartGrid = 'OFF';
    Mode = 'New';%_Sterile';
    Twin_mNu4Sq   = 0;%1e2;
    Twin_sin2T4   = 0;%0.06;
    
    extraStr = '';
    switch Mode
        case 'Old'
            RandMC = [1:151,500:643]*1e3;
            SysBudget =22;
        case 'New'
            SysBudget =24;
            if range== 95
                RandMC = 1:1500;
            elseif range==65
                RandMC = [1:654,868:931,1219:1500];
            elseif range==40
                RandMC = 1:2000;
            end
        case 'New_Sterile'
            SysBudget =24;
            RandMC = 1:1000;
            extraStr = sprintf('_m4Sq%.0feV2_sin2t4%.2f',Twin_mNu4Sq,Twin_sin2T4);
    end
    
    %% set up model for contour calculation 
  
%% configure RunAnalysis object
Real = MultiRunAnalysis('RunList','KNM1',...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'fixPar',freePar,...
    'NonPoissonScaleFactor',1,...
    'SysBudget',200,...
    'minuitOpt','min ; minos',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AngularTFFlag','ON',...
    'SynchrotronFlag','ON',...
    'RadiativeFlag','ON',...
    'DopplerEffectFlag','FSD',...
    'BKG_PtSlope',-2.2*1e-06);
Real.exclDataStart = Real.GetexclDataStart(range);
%% configure Sterile analysis object
SterileArg = {'RunAnaObj',Real,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};

S = SterileAnalysis(SterileArg{:});

    %% init
    nContours = numel(RandMC);
    mnu4Sq   = cell(nContours,1);
    sin2T4   = cell(nContours,1);
    chi2     = cell(nContours,1);
    chi2_ref = cell(nContours,1);
    savefile = cell(nContours,1);
    FitResults_ref = cell(nContours,1);
    DeltaChi2 = zeros(nContours,1);
    
      mNu4Sq_contour = cell(nContours,1);
    sin2T4_contour = cell(nContours,1);
    %% load grid (or calculate if doesn't exist)
    mnu4Sq_bf = zeros(numel(RandMC),1);
    sin2T4_bf = zeros(numel(RandMC),1);
    chi2min_bf   = zeros(numel(RandMC),1);
    chi2min_null  = zeros(numel(RandMC),1);
    
    for i=RandMC
        progressbar(i/nContours)
        [mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile{i},FitResults_Null] = KSN1GridSearch('range',range,...
            'nGridSteps',nGridSteps,...
            'chi2',chi2Str,...
            'DataType',DataType,...
            'freePar','E0 Bkg Norm',...
            'RunList',RunList,...
            'SmartGrid',SmartGrid,...
            'RecomputeFlag','OFF',...
            'RandMC',i,...
            'SysBudget',SysBudget,...
            'Twin_mNu4Sq',Twin_mNu4Sq,...
            'Twin_sin2T4',Twin_sin2T4);
        
        chi2tmp   = chi2{i};
        
        % find best fit
        mnu4Sqtmp = mnu4Sq{i};
        sin2T4tmp = sin2T4{i};
        
%         [row, col] = find(chi2tmp == min(chi2tmp(:)));
%         mnu4Sq_bf(i) =  mnu4Sqtmp(col,row);
%         sin2T4_bf(i)  = sin2T4tmp(col,row);
%         chi2min_bf(i) = min(chi2tmp(:));
         chi2min_null(i) = FitResults_Null.chi2min;
%         DeltaChi2(i) = chi2min_null(i)-min(min(chi2tmp));
        
        S.InterpMode = InterpMode;
        S.sin2T4 = sin2T4tmp;
        S.chi2 = chi2tmp;
        S.mNu4Sq = mnu4Sqtmp;
        S.mNuSq = zeros(25,25);
        S.E0 = zeros(25,25);
        S.Interp1Grid;
        S.ContourPlot('BestFit','ON'); close;
        
        mNu4Sq_contour{i} = S.mNu4Sq_contour;
        sin2T4_contour{i} = S.sin2T4_contour;
        mnu4Sq_bf(i) = S.mNu4Sq_bf;
        sin2T4_bf(i) = S.sin2T4_bf;
        chi2min_bf(i) = S.chi2_bf;
       
        DeltaChi2(i) = chi2min_null(i)-chi2min_bf(i);
    end
    %
    mnu4Sq_bf = mnu4Sq_bf(RandMC);
    sin2T4_bf =sin2T4_bf(RandMC);
    chi2min_bf = chi2min_bf(RandMC);
    chi2min_null = chi2min_null(RandMC);
    DeltaChi2 = DeltaChi2(RandMC);
    
    d = importdata(savefile{200});
    dof = d.FitResults_Null.dof-2;
    %% significant best fits
    
    ClosedLog95 =  DeltaChi2>= GetDeltaChi2(0.95,2);
    ClosedFrac95 = sum(ClosedLog95)/nContours;
    RelErr = @(n,p) sqrt((n*p*(1-p)))/n;
    
    fprintf('%.0f%% C.L. : fraction of significant best fits = %.1f  +- %.1f  (%.0f out of %.0f)\n',...
        95,ClosedFrac95*100,RelErr(nContours,ClosedFrac95)*100,sum(ClosedLog95),nContours);
    
    save(savename,'RandMC','mnu4Sq_bf','sin2T4_bf','chi2min_bf','chi2min_null','DeltaChi2','dof',...
        'ClosedLog95','ClosedFrac95','RelErr','mNu4Sq_contour','sin2T4_contour')
end
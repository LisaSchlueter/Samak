
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
% 
% Configuration for First Tritium May

% P. Morales 2018
% Last update 23/04/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%Initialization
clear; close all;
addpath(genpath('../../../Samak2.0'));

%Parameters
nfit            = 1;
runtime        = (1)*60*60*24;
%runtime        = 210.2189781*18;

nTeBinFac       = 5;
TD              = 'Run40531';
Fitter          = 'matlab';
FPD_Pixel       = 1:148;

FSD             = 'DOSS';
BCK             = 'XmasData';
Mode            = 'Read';

% Computing
savespectrum    = 'OFF';
FitFlag         = 'OFF';
UseCovMat       = 'OFF';

% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Initialization
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute'}; %MACE+WGTSIS %SSCW_DCOfficial %SSCW_DCperqU

opt_bkg = {...
    'BKG_Flag',BCK};

opt_fsd= {'DTFSD',FSD};

opt_katrin = {...
    'TD',TD,...
    'TimeSec',runtime,...
    'Mode',Mode};

opt_fpd = {...
    'FPD_Pixel',FPD_Pixel};

A = ref_ftmultipix(...
    opt_katrin{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_bin{:},...
    opt_fpd{:},...
    'mnuSq_i',mnuSq_t,...
    'UseParallelRF','OFF',...
    'FPD_Segmentation','MULTIPIXEL');

if strcmp(UseCovMat,'ON')
    if exist('results/CM.mat','file') == 2
        CM = load('results/CM.mat');
        MegaCovMatFrac = CM.MegaCovMatFrac;
    else
        A_CM = ref_ftmultipix(...
            opt_katrin{:},...
            opt_wgtsmace{:},...
            opt_fsd{:},...
            opt_bkg{:},...
            opt_bin{:},...
            'mnuSq_i',mnuSq_t,...
            'UseParallelRF','OFF',...
            'FPD_Segmentation','OFF');
        
        %2. Define your CM Class with your number of trials, choose SysEffect,...
        SysEffect = struct(...
            'RF_EL','OFF',...  % Response Function(RF) EnergyLoss
            'RF_BF','OFF',...  % R F B-Fields
            'RF_RX','OFF',...  % RF Column Density, Cross Section
            'FSD','ON',...
            'TASR','ON');     % Tritium Activity Stack Runs
        C = CovarianceMatrix('StudyObject',A_CM,'nTrials',1500,...
            'RecomputeFlag','ON',...
            'RunFlag','OFF',...
            'SysEffect',SysEffect,...
            'MACE_Ba_T_RelErr',1e-2,...
            'WGTS_B_T_RelErr',1e-2,...
            'WGTS_CD_MolPerCm2_RelErr',0.1,...
            'ISXsection_RelErr',1e-2);
        %C.WGTS_TASR_RelErr = 0.01;
        
        %%4.Compute Covariance Matrix for your SysEffect; here RF (Response Function)
        C.ComputeCM_FSD;

        
        % Build mega Covariance Matrix
        CovMatpixCell = repmat({C.CovMatFrac},1,length(1));
        MegaCovMatFrac = blkdiag(CovMatpixCell{:});
    end
    
    A.ComputeTBDDS();
    A.ComputeTBDIS();
    
    MegaCovMat = A.TBDIS(:).*MegaCovMatFrac.*A.TBDIS(:)';
    
    % Define Covariance MatriX for fit
    FitCovErrMatrix = diag(A.TBDIS(:)) + MegaCovMat;
    
end
%%
    A.ComputeTBDDS();
    A.ComputeTBDIS();
    A.AddStatFluctTBDIS();
    Data = [...
        A.qU, ...
        A.TBDIS, ...
        A.TBDISE];
    Afit = ref_ftmultipix(...
        opt_katrin{:},...
        opt_wgtsmace{:},...
        opt_fsd{:},...
        opt_fpd{:},...
        'BKG_Flag','ON',...
        'BKG_Type','FLAT',...
        'BKG_RateAllFPDSec',0,...
        opt_bin{:},...
        'mnuSq_i',mnuSq_t,...
        'UseParallelRF','ON',...
        'FPD_Segmentation','OFF',...
        'FPD_Pixel',FPD_Pixel);
    Afit.ComputeTBDDS();
    Afit.ComputeTBDIS();  
    
    F = FITC('SO',Afit,'DATA',Data,'fitter',Fitter,'pulls',[4,inf(1,Afit.nPixels*2+1)],'savepdf','OFF');
    
    
    
    
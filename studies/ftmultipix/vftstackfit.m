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

DS = load('stackDTallmpix');

TD = 'stackRunDTallmpix';
%Parameters
nfit            = 1;

nTeBinFac       = 5;
Fitter          = 'matlab';
FPD_Pixel       = 1;

FSD             = 'DOSS';
BCK             = 'FLAT';
Mode            = 'Read';
Segment         = 'SINGLEPIXEL';

% Computing
FitFlag         = 'ON';
UseCovMat       = 'ON';

% Parametrization: True Value
mnuSq_t         = (0.00)^2;

% Initialization
opt_bin = {...
    'nTeBinningFactor',nTeBinFac};

opt_wgtsmace = {...
    'KTFFlag','Compute'};

opt_bkg = {...
    'BKG_Flag',BCK};

opt_fsd= {'DTFSD',FSD,...
    'WGTS_MolFrac_DT',DS.WGTS_MolFrac_DT};

opt_katrin = {...
    'TD',TD,...
    'TimeSec',DS.TimeSec,...
    'Mode',Mode,...
    'WGTS_CD_MolPerCm2',DS.WGTS_CD_MolPerCm2};



if strcmp(UseCovMat,'ON')
    if exist('results/VFT_1500Trials.mat','file') == 2
        CM = load('results/CombiCMVFT.mat');
        CovMatFrac = CM.CMcombifrac(1:24,1:24);
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
        
        %4.Compute Covariance Matrix for your SysEffect; here RF (Response Function)
        C.ComputeCM_FSD;

    end

    
end
%%

% pixlist = [1;2;3;7;11;12;14;15;16;17;20;21;22;24;25;28;29;30;31;32;33;36;37;41;...
%     42;43;46;50;53;55;56;60;61;62;63;64;66;67;68;69;71;73;74;77;78;79;81;...
%     82;87;88;89;91;92;93;94;98;100;101;102;103;104;107;108;109;111;113;...
%     114;115;116;117;118;119;120;124;126;127;128;131;132;134;136;146];

% without outliers and without outer two rings
% pixlist = [1;2;3;7;11;12;14;15;16;17;20;21;22;24;25;28;31;32;33;36;37;41;...
%     42;43;46;50;53;55;56;60;61;62;63;64;66;67;68;69;71;73;74;77;78;79;81;...
%     82;87;88;89;91;92;93;94;98;100;101;102;103;104;107;108;109;111;113;...
%     114;115;116;117;118;119;120];

%pixlist = ([1:112,114:123])';
pixlist = (1:5)';
%for p = 1:length(pixlist)
% pixlist = (1:76)';
% par = zeros(148,4); 
% err = zeros(148,4); 
% chi2min = zeros(148,1); 
% CMfit = cell(148,1);
for p = 1:length(pixlist)
    
    pix = pixlist(p);
   
    opt_fpd = {...
    'pixel',pix,...
    'FPD_Segmentation','SINGLEPIXEL'};

Data = [DS.qU(:,pix),DS.TBDIS(:,pix),sqrt(DS.TBDIS(:,pix))];

Afit = ref_stack(...
    opt_katrin{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    'BKG_Flag','ON',...
    'BKG_Type',BCK,...
    'BKG_RateAllFPDSec',0,...
    opt_bin{:},...
    'mnuSq_i',mnuSq_t,...
    opt_fpd{:});

Afit.ComputeTBDDS();
Afit.ComputeTBDIS();

% Build mega Covariance Matrix
%CovMatpixCell = repmat({CovMatFrac},1,Afit.nPixels);
%MegaCovMatFrac = blkdiag(CovMatpixCell{:}) + diag(1./Data(:,2));


if strcmp(FitFlag,'ON')
    F = FITC('SO',Afit,'DATA',Data,'fitter',Fitter,'pulls',[4,inf(1,Afit.nPixels*2+1)],...
        'savepdf','OFF','exclDataStart',1,'exclDataStop',24,'chi2name','chi2P');%,...
      %  'COVMATFrac',MegaCovMatFrac);
    par(pix,:) = F.RESULTS{1}; err(pix,:) = F.RESULTS{2}; chi2min(pix,:) = F.RESULTS{3};
    CMfit{pix}= F.RESULTS{5};
end

end % for pixels

TOTALRESULTS = {par,err,chi2min};
F.RESULTS = TOTALRESULTS;
%F.savepdffunction()
F.savepdfiterative(pixlist,DS.qU(:,pixlist),DS.TBDIS(:,pixlist),sqrt(DS.TBDIS(:,pixlist)),CMfit)

   
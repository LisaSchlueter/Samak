% --------------------------------------------------------------------------------
% Goal: 
% check if response function broadening can compensate for different qU when stacking runs
% using Twin MC
% - all runs have common endpoint
% - other slow control identical to data (subrun wise qU, column density,T2,...)
%
% Method:
% - calculate twins
% - calculate model with broadened RF
% - fit with only 1 common endpoint
% - estimate neutrino mass bias
% 
% Lisa Schlüter
% January 2020
% --------------------------------------------------------------------------------
%%  common settings
RecomputeFlag = 'ON';
RunList = 'KNM2_RW1';    % all runs
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Twin';     
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel';     % uniform FPD
range = 40;                 % in eV below E0
chi2 = 'chi2Stat';
TwinBias_Q = 18573.70;
Twin_SameqUFlag = 'OFF';

%% load if possible
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
savename = [savedir,sprintf('knm2_RunStackingImp_%s_%.0feV.mat',RunList,range)];

if strcmp(Twin_SameqUFlag,'ON') 
    savename = strrep(savename,'.mat','_SameqU.mat');
end

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    %  collect common arguments                               
    CommonArg = {'RunList',RunList,...
        'fixPar',fixPar,...
        'DataType',DataType,...
        'FSDFlag',FSDFlag,...
        'ELossFlag',ELossFlag,...
        'AnaFlag',AnaFlag,'chi2',chi2,...
        'TwinBias_Q',TwinBias_Q,...
        'Twin_SameqUFlag',Twin_SameqUFlag};
   
    %% 1. Twins with same endpoint
    A = MultiRunAnalysis(CommonArg{:});  % constant endpoint
    A.exclDataStart = A.GetexclDataStart(range);              % get correct range
    A.Fit;                                                    % fit
    FitResult = A.FitResult;                                  % save

    RFsigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = RFsigma;
    A.ModelObj.InitializeRF('RebinMode','Integral'); 
    A.Fit;
    FitResultImp = A.FitResult;   
    
    MakeDir(savedir); %create directory if it necessary
    save(savename,'FitResult','FitResultImp','CommonArg','A');
end

%% display neutrino mass shift
fprintf('Mean');
fprintf('------------------------------------------- \n')
fprintf('mNuSq bias =  %.3f eV^2 (reg.) vs. %.3f eV^2 (imp.)  \n',FitResult.par(1),FitResultImp.par(1));
fprintf('E0    bias =  %.3f eV vs. %.3f eV (imp.) \n',FitResult.par(2)+A.ModelObj.Q_i-TwinBias_Q,FitResultImp.par(2)+A.ModelObj.Q_i-TwinBias_Q);
fprintf('------------------------------------------- \n')


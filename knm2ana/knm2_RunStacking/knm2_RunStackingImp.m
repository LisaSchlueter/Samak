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
RunList = 'KNM2_Prompt';    % all runs
range = 40;                 % in eV below E0
TwinBias_Q = 18573.7;

%% load if possible
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
savename = [savedir,sprintf('knm2_RunStackingImp_%s_%.0feV_E0%.2feV.mat',RunList,range,TwinBias_Q)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    %  collect common arguments                               
    CommonArg = {'RunList',RunList,...
        'fixPar','mNu E0 Bkg Norm',...
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'TwinBias_Q',TwinBias_Q,...
        'NonPoissonScaleFactor',1};
    
    % twin runs all have the same MTD
    B = MultiRunAnalysis(CommonArg{:},'Twin_SameqUFlag','ON');   
    B.exclDataStart = B.GetexclDataStart(range);                           
    B.Fit;
    FitResult_sameqU = B.FitResult;
    
    % twin runs have MTDs from real data
    A = MultiRunAnalysis(CommonArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    A.Fit;       % fit without RF broadening
    FitResult = A.FitResult;
    
    RFsigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = RFsigma;
    A.ModelObj.InitializeRF('RebinMode','Integral');
    A.Fit;     % fit with RF broadening
    FitResultImp = A.FitResult;
    
    MakeDir(savedir); %create directory if it necessary
    save(savename,'FitResult','FitResultImp','FitResult_sameqU','CommonArg','A');
end

%% display neutrino mass shift
fprintf('------------------------------------------- \n')
fprintf('------------- neutrin mass bias ---------- \n')
fprintf('regular twins     : %.4f eV^2  \n',FitResult.par(1));
fprintf('regular twins imp : %.4f eV^2  \n',FitResultImp.par(1));
fprintf('qU same twins     : %.4f eV^2  \n',FitResult_sameqU.par(1));
fprintf('------------- Endpoint bias --------------- \n')
fprintf('regular twins     : %.3f eV  \n',FitResult.par(2)+A.ModelObj.Q_i-TwinBias_Q);
fprintf('regular twins imp : %.3f eV  \n',FitResultImp.par(2)+A.ModelObj.Q_i-TwinBias_Q);
fprintf('qU same twins     : %.3f eV  \n',FitResult_sameqU.par(2)+A.ModelObj.Q_i-TwinBias_Q);
fprintf('------------------------------------------- \n')


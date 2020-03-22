% Fit KNM2 alternative runlist

freePar = 'E0 Bkg Norm';
DataType = 'Real';
RunList = 'KNM2_Prompt';
range = 40;                % fit range in eV below endpoint
AltPixList = 'AziHalfNS';%'AziHalfEW';  % defines alternative pixel list
RecomputeFlag = 'OFF';

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_PixListAlt_%s_%s_%s_%.0feV.mat',...
    savedir,AltPixList,DataType,strrep(freePar,' ',''),range);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
RunAnaArg = {'RunList',RunList,...     % define run number -> see GetRunList
    'fixPar',freePar,...               % free Parameter !!
    'DataType',DataType,...            % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.7,...
    'ROIFlag','14keV',...
    'DopplerEffectFlag','FSD',...
    'RingMerge',AltPixList};

M = MultiRunAnalysis(RunAnaArg{:});          % init model, read data
M.exclDataStart = M.GetexclDataStart(range); % set fit range

R = RingAnalysis('RunAnaObj',M,'RingList',M.RingList);
R.FitRings;

FitResult = R.FitResult(1);
PixList = arrayfun(@(x) x.PixList,R.MultiObj,'UniformOutput',false);
RunList = M.RunList;
Q_i = M.ModelObj.Q_i;

mNuSq    = FitResult.par(:,1); 
mNuSqErr = FitResult.err(:,1);
E0       = FitResult.par(:,2)+Q_i;
E0Err    = FitResult.err(:,2);

MakeDir(savedir);
save(savename,'FitResult','PixList','RunList','Q_i','mNuSq','E0','mNuSqErr','E0Err');
end
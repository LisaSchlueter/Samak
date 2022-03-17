% Fit KNM2 alternative runlist

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
AltRunList = 'KNM2_Prompt';       % defines alternative pixel list
RecomputeFlag = 'ON';
FSDFlag = 'KNM2_0p1eV';
BKG_PtSlope = 3*1e-06;
% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_AltRunList_%s_%s_%s_%.0feV_%s_BkgPt%.2g.mat',...
    savedir,AltRunList,DataType,strrep(freePar,' ',''),range,FSDFlag,BKG_PtSlope*1e6);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
    fprintf('load result from %s \n',savename);
else
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList',AltRunList,...     % define run number -> see GetRunList
        'fixPar',freePar,...               % free Parameter !!
        'DataType',DataType,...            % Real, Twin or Fake
        'FSDFlag',FSDFlag,...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1.112,...
        'TwinBias_Q',18573.7,...
        'DopplerEffectFlag','FSD',...
        'RingMerge','Full',...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'BKG_PtSlope',BKG_PtSlope};
    
    M = MultiRunAnalysis(RunAnaArg{:});          % init model, read data
    M.exclDataStart = M.GetexclDataStart(range); % set fit range
    
    M.Fit;
    FitResult = M.FitResult;
    
    RunList = M.RunList;
    Q_i = M.ModelObj.Q_i;
    
    mNuSq    = FitResult.par(1);
    mNuSqErr = FitResult.err(1);
    E0       = FitResult.par(2)+Q_i;
    E0Err    = FitResult.err(2);
    
    MakeDir(savedir);
    save(savename,'FitResult','RunList','Q_i','mNuSq','E0','mNuSqErr','E0Err','RunAnaArg');
    fprintf('save result to %s \n',savename);
end
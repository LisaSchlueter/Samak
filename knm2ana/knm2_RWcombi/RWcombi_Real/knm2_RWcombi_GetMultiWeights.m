function MultiWeights = knm2_RWcombi_GetMultiWeights

RecomputeFlag = 'OFF';
savedir = [getenv('SamakPath'),'inputs/Plasma/Knm2/'];
savename = [savedir,'Knm2_PlasmaRWcombi_MultiWeights.mat'];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,'MultiWeights');
else
    % Get time spent in each RW period
    RunAnaArg = {'fixPar','E0 Bkg Norm',...         % free Parameter !!
        'DataType','Real',...              % Real, Twin or Fake
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'NonPoissonScaleFactor',1};
    
    % build object of MultiRunAnalysis class
    A1 = MultiRunAnalysis('RunList','KNM2_RW1',RunAnaArg{:});
    A2 = MultiRunAnalysis('RunList','KNM2_RW2',RunAnaArg{:});
    A3 = MultiRunAnalysis('RunList','KNM2_RW3',RunAnaArg{:});
    
    MultiWeights = zeros(3,1);
    MultiWeights(1) = A1.RunData.TimeSec;
    MultiWeights(2) = A2.RunData.TimeSec;
    MultiWeights(3) = A3.RunData.TimeSec;
    MultiWeights = MultiWeights./sum(MultiWeights);
    
    %save
    MakeDir(savedir);
    save(savename,'MultiWeights','RunAnaArg');
end
end
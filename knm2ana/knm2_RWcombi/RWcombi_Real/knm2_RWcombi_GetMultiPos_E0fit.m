function MultiPos = knm2_RWcombi_GetMultiPos_E0fit
% get Multi pos from E0 fit (preliminary method)
RecomputeFlag = 'OFF';
savedir = [getenv('SamakPath'),'inputs/Plasma/Knm2/'];
savename = [savedir,'Knm2_PlasmaRWcombi_MultiPos_E0fit.mat'];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,'MultiPos');
else
    % Get time spent in each RW period
    RunAnaArg = {'fixPar','E0 Bkg Norm',...         % free Parameter !!
        'DataType','Real',...              % Real, Twin or Fake
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'NonPoissonScaleFactor',1,...
        'exclDataStart',1}; % 90 eV range
    
    % build object of MultiRunAnalysis class
    A1 = MultiRunAnalysis('RunList','KNM2_RW1',RunAnaArg{:});
    A2 = MultiRunAnalysis('RunList','KNM2_RW2',RunAnaArg{:});
    A3 = MultiRunAnalysis('RunList','KNM2_RW3',RunAnaArg{:});
    
    A1.Fit;
    A2.Fit;
    A3.Fit;
    
    MultiPos = zeros(3,1);
    MultiPos(1) = A1.FitResult.par(2);
    MultiPos(2) = A2.FitResult.par(2);
    MultiPos(3) = A3.FitResult.par(2);
    MultiPos = MultiPos-mean(MultiPos);
    
    %save
    MakeDir(savedir);
    save(savename,'MultiPos','RunAnaArg');
end
end
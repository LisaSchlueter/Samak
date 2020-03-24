% Uniform fit on KNM2 data
% random half of all golden runs
% March 2020, Lisa

nFits = 500;
freePar = 'E0 Bkg Norm';
DataType = 'Real';
range = 40;               % fit range in eV below endpoint
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_RunListRandHalf_%s_%s_%.0feV_%.0ffits.mat',...
            savedir,DataType,strrep(freePar,' ',''),range,nFits);


if exist(savename,'file')
    load(savename,'RunList','FitResult','RunAnaArg')
else
    
    RunList = cell(nFits,1);
    FitResult = cell(nFits,1);
    RunAnaArg = {'RunList','KNM2_RandHalf',...  % define run number -> see GetRunList
        'fixPar',freePar,...         % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'ROIFlag','14keV',...
        'DopplerEffectFlag','FSD'};%,...
    %'Twin_SameqUFlag','ON'};
    
    for i=1:nFits
        %% build object of MultiRunAnalysis class
        D = MultiRunAnalysis(RunAnaArg{:});
        
        % modify some parameters in your analysis
        D.exclDataStart = D.GetexclDataStart(range); % find correct data, where to cut spectrum
        
        %% Fit -> fit results are in property: A.FitResult
        D.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        D.Fit;
        
        % save
        RunList{i} = D.RunList;
        FitResult{i} = D.FitResult;
    end
    MakeDir(savedir);
    save(savename,'RunList','FitResult','RunAnaArg');
end

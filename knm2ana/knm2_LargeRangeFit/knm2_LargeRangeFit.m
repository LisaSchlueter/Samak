
function [E0, E0Err] = knm2_LargeRangeFit
% uniform fit to large range
% get rear wall period-wise endpoint
% [-90, -45] --> only subruns, that are not used for nu-mass analysis
RecomputeFlag = 'OFF';
chi2 = 'chi2Stat';
%% load or calc
savedir = [getenv('SamakPath'),'knm2ana/knm2_LargeRangeFit/results/'];
savename = sprintf('%sknm2FS_LargeRangeFit_%s.mat',savedir,chi2);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
    fprintf('Load large range fits from file %s \n',savename);
else
    RunAnaArg = {...
        'fixPar','E0 Bkg Norm',...         % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2',chi2,...                    % statistics only
        'ROIFlag','14keV',...
        'SysBudget',33,...
        'exclDataStart',1,...
        'exclDataStop',10};
    if strcmp(chi2,'chi2Stat')
        RunAnaArg = [RunAnaArg,'NonPoissonScaleFactor',1];
    else
        RunAnaArg = [RunAnaArg,'NonPoissonScaleFactor',1.112];
    end
    %% build object of MultiRunAnalysis class
    nPeriods = 3;
    nRuns = [171,97,93];
    E0mean = 18573.69;% get back twins with (overall mean of 361 runs)
    FitResults = cell(nPeriods,1);
    for i=1:nPeriods
        clear A;
        A = MultiRunAnalysis('RunList',sprintf('KNM2_RW%.0f',i),RunAnaArg{:},...
            'TwinBias_Q',E0mean*ones(nRuns(i),1));
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
        A.Fit;
        FitResults{i} = A.FitResult;
    end
    Q_i = A.ModelObj.Q_i;
    MakeDir(savedir);
    %A.PlotFit('FitResultsFlag','OFF','MaxBkgRange',-40);
    save(savename,'FitResults','E0mean','nRuns','RunAnaArg','nPeriods','Q_i');
    fprintf('Save large range fits from file %s \n',savename);
end

%%
E0    = cell2mat(cellfun(@(x) x.par(2)+Q_i,FitResults,'UniformOutput',0));
E0Err = cell2mat(cellfun(@(x) x.err(2),FitResults,'UniformOutput',0));

E0twins = knm2FS_GetE0Twins('SanityPlot','OFF');
fprintf('E0 std from twins : %.0f meV \n',std(E0twins)*1e3);
fprintf('E0 std from lr fit: %.0f meV \n',std(E0)*1e3);
end

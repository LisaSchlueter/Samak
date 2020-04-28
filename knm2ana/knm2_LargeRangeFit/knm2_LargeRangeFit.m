
function [E0, E0Err] = knm2_LargeRangeFit(varargin)

% uniform fit to large range
% get rear wall period-wise endpoint
% [-90, -45] --> only subruns, that are not used for nu-mass analysis
% one additional background point at +135 eV
p = inputParser;
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat','chi2CMShape'}));
p.addParameter('DataType','Real',@(x)ismember(x,{'Twin','Real'}));
p.addParameter('NonPoissonScaleFactor',1.112,@(x)isfloat(x));
p.parse(varargin{:});

RecomputeFlag = p.Results.RecomputeFlag;
chi2          = p.Results.chi2;
DataType      = p.Results.DataType;
NonPoissonScaleFactor = p.Results.NonPoissonScaleFactor;

%% load or calc
savedir = [getenv('SamakPath'),'knm2ana/knm2_LargeRangeFit/results/'];
savename = sprintf('%sknm2FS_LargeRangeFit_%s_%s_NPfactor%.3f.mat',savedir,DataType,chi2,NonPoissonScaleFactor);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
    fprintf('Load large range fits from file %s \n',savename);
else
    RunAnaArg = {...
        'fixPar','E0 Bkg Norm',...         % free Parameter !!
        'DataType',DataType,...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2',chi2,...                    % statistics only
        'ROIFlag','Default',...
        'SysBudget',35,...
        'exclDataStart',[1:10,38],...   
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor};

    %% build object of MultiRunAnalysis class
    nPeriods = 3;
    nRuns = [171,97,93];
    FitResults = cell(nPeriods,1);
    for i=1:nPeriods
        clear A;
        A = MultiRunAnalysis('RunList',sprintf('KNM2_RW%.0f',i),RunAnaArg{:});
        A.Fit;
        FitResults{i} = A.FitResult;
    end
    Q_i = A.ModelObj.Q_i;
    MakeDir(savedir);
   % A.PlotFit('FitResultsFlag','OFF','MaxBkgRange',136);
    save(savename,'FitResults','nRuns','RunAnaArg','nPeriods','Q_i');
    fprintf('Save large range fits from file %s \n',savename);
end

%%
E0    = cell2mat(cellfun(@(x) x.par(2)+Q_i,FitResults,'UniformOutput',0));
E0Err = cell2mat(cellfun(@(x) x.err(2),FitResults,'UniformOutput',0));

fprintf('E0 std from lr fit: %.0f meV \n',std(E0)*1e3);
end

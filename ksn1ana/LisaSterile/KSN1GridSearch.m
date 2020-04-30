function [mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch(varargin)
% parallel grid search for sterile ksn1 analysis
% Lisa, April 2020

p = inputParser;
p.addParameter('range',40,@(x)isfloat(x));
p.addParameter('nGridSteps',50,@(x)isfloat(x)); % time on server: 25 points ~ 7-10 minutes, 50 points ~ 40 minutes on csltr server
p.addParameter('chi2','chi2CMShape',@(x)ismember(x,{'chi2CMShape','chi2Stat'}));
p.addParameter('DataType','Twin',@(x)ismember(x,{'Twin','Real'}));
p.addParameter('freePar','E0 Bkg Norm',@(x)ischar(x)); % check ConvertFixPar.m for options
p.addParameter('RunList','KNM1',@(x)ischar(x)); 
p.addParameter('SmartGrid','OFF',@(x)ismember(x,{'ON','OFF'})); % work in progress
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'})); 

p.parse(varargin{:});

range         = p.Results.range;
nGridSteps    = p.Results.nGridSteps;
chi2          = p.Results.chi2;
DataType      = p.Results.DataType;
freePar       = p.Results.freePar;
RunList       = p.Results.RunList;
SmartGrid     = p.Results.SmartGrid;
RecomputeFlag = p.Results.RecomputeFlag;

if strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor=1.064;
else
    NonPoissonScaleFactor=1;
end

%% label
if strcmp(SmartGrid,'ON')
    AddSin2T4 = 0.1;
    extraStr = sprintf('_SmartGrid%.0e',AddSin2T4);
else
    extraStr = '';
end

savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
MakeDir(savedir);

savefile = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
    savedir,RunList,DataType,strrep(freePar,' ',''),range,chi2,nGridSteps,extraStr);
%% load or calculate grid
if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile,'mnu4Sq','sin2T4','chi2','chi2_ref')
    fprintf('load grid from file %s \n',savefile)
else
    if range<=40
        FSDFlag = 'Sibille0p5eV';
    elseif range>40
        FSDFlag = 'SibilleFull';
    end
    
    RunAnaArg = {'RunList',RunList,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'TwinBias_Q',18573.73};
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    if strcmp(DataType,'Real')
        T.fixPar = [freePar,'mnu4Sq , sin2T4'];
        T.InitFitPar;
        T.pullFlag = 9; % limit sin4 to [0,0.5] and m4 [0 range^2]
    end
    
    T.Fit;
    FitResults_ref = T.FitResult;
    chi2_ref       = T.FitResult.chi2min;
    
    if strcmp(DataType,'Real')
        T.fixPar = freePar;
        T.InitFitPar;
        T.pullFlag = 99; % means no pull
    end
    %% define msq4 - sin2t4 grid
    switch SmartGrid
        case 'OFF'
            sin2T4      = logspace(-3,log10(0.5),nGridSteps); %linspace(0.001,0.5,nGridSteps)
            mnu4Sq      = logspace(0,log10((range+5)^2),nGridSteps)';
            mnu4Sq      = repmat(mnu4Sq,1,nGridSteps);
            sin2T4      = repmat(sin2T4,nGridSteps,1);
        case 'ON'
            [mnu4Sq,sin2T4] = GetSmartKsn1Grid('range',range,'ConfLevel',95,...
                'nGridSteps',nGridSteps,'AddSin2T4',AddSin2T4,...
                'SanityPlot','ON');
    end
    %% make copy of model for parallel computing
    D = copy(repmat(T,nGridSteps.*nGridSteps,1));
    
    %% scan over msq4-sin2t4 grid
    D              = reshape(D,numel(D),1);
    chi2Grid       = zeros(nGridSteps*nGridSteps,1);
    FitResultsGrid = cell(nGridSteps*nGridSteps,1);
    mnu4Sq_Grid    = reshape(mnu4Sq',nGridSteps*nGridSteps,1);
    sin2T4_Grid    = reshape(sin2T4',nGridSteps*nGridSteps,1);
    
    parfor i= 1:(nGridSteps*nGridSteps)
        D(i).SimulateStackRuns;
        D(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
        D(i).Fit
        chi2Grid(i) = D(i).FitResult.chi2min;
        FitResultsGrid{i} = D(i).FitResult;
    end
    
    chi2       = reshape(chi2Grid,nGridSteps,nGridSteps);
    FitResults = reshape(FitResultsGrid,nGridSteps,nGridSteps);
    %% save
    save(savefile,'chi2_ref','FitResults_ref','RunAnaArg',...
        'chi2','mnu4Sq','sin2T4','FitResults');
    fprintf('save file to %s \n',savefile)
    
    try 
        %% get contour at X sigma (can also be calculated later for any CL)
        [mnu4Sq_contour95, sin2T4_contour95] = KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,95);
        save(savefile,'mnu4Sq_contour95','sin2T4_contour95','-append')
        [mnu4Sq_contour90, sin2T4_contour90] = KSN1Grid2Contourr(mnu4Sq,sin2T4,chi2,chi2_ref,90);
        save(savefile,'mnu4Sq_contour90','sin2T4_contour90','-append')
    catch
    end
end
end


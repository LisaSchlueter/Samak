function [mnu4Sq,sin2T4,chi2,chi2_ref,savefile,FitResults_Null] = KSNXGridSearch(varargin)
% parallel grid search for sterile ksnx analysis
% Lisa, April 2020

p = inputParser;
p.addParameter('range',40,@(x)isfloat(x));
p.addParameter('nGridSteps',50,@(x)isfloat(x)); % time on server: 25 points ~ 7-10 minutes, 50 points ~ 40 minutes on csltr server
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2CMShape','chi2Stat'}));
p.addParameter('DataType','Fake',@(x)ismember(x,{'Twin','Real','Fake'}));
p.addParameter('freePar','E0 Bkg Norm',@(x)ischar(x)); % check ConvertFixPar.m for options
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SysEffect','all',@(x)ischar(x)); % if chi2CMShape: all or only 1
p.addParameter('RandMC','OFF',@(x)ischar(x) || isfloat(x)); % randomize twins if RandMC is float
p.addParameter('pullFlag',99,@(x)isfloat(x)); % if above 12 --> no pull
p.addParameter('SysBudget',66,@(x)isfloat(x));
p.addParameter('ELossFlag','KatrinT2A20',@(x)ischar(x));
p.addParameter('AngularTFFlag','ON',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});

range         = p.Results.range;
nGridSteps    = p.Results.nGridSteps;
chi2          = p.Results.chi2;
DataType      = p.Results.DataType;
freePar       = p.Results.freePar;
RecomputeFlag = p.Results.RecomputeFlag;
SysEffect     = p.Results.SysEffect;
RandMC        = p.Results.RandMC;
pullFlag      = p.Results.pullFlag;
SysBudget     = p.Results.SysBudget;
ELossFlag     = p.Results.ELossFlag;
AngularTFFlag = p.Results.AngularTFFlag;

if strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor=1;
else
    NonPoissonScaleFactor=1;
end

FakeInitFile = @ref_KSNX_KATRIN_Final;

%% label
extraStr = '';
if strcmp(chi2,'chi2CMShape')
    extraStr = [extraStr,sprintf('_Budget%.0f',SysBudget)];
    if ~strcmp(SysEffect,'all')
        extraStr = [extraStr,sprintf('_%s',SysEffect)];
    end
end

if ~strcmp(ELossFlag,'KatrinT2')
    extraStr = [extraStr,sprintf('_%s',ELossFlag)];
end

if ~strcmp(AngularTFFlag,'OFF')
    extraStr = [extraStr,'_AngTF'];
end

if isfloat(RandMC) && strcmp(DataType,'Twin')
    extraStr = sprintf('%s_RandMC%.0f',extraStr,RandMC);
    savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/RandomizedMC/'];
else
    savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
end

if pullFlag<=12
    extraStr = sprintf('%s_pull%.0f',extraStr,pullFlag);
end
MakeDir(savedir);

savefile = sprintf('%sKSNX_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
    savedir,func2str(FakeInitFile),DataType,strrep(freePar,' ',''),range,chi2,nGridSteps,extraStr);

%% check if a file with larger n grid size exist. If yes, load that on instead
if ~exist(savefile,'file') && strcmp(RecomputeFlag,'OFF') && nGridSteps<50
    nGridSteps_i = nGridSteps;
    nGridSteps = 50;
    savefile = sprintf('%sKSNX_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
        savedir,func2str(FakeInitFile),DataType,strrep(freePar,' ',''),range,chi2,nGridSteps,extraStr);
    nGridSteps = nGridSteps_i;
    if ~exist(savefile,'file') % if grid for 50 also doesn't exist -> go back to 25 grid
        savefile = sprintf('%sKSNX_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid%s.mat',...
            savedir,func2str(FakeInitFile),DataType,strrep(freePar,' ',''),range,chi2,nGridSteps,extraStr);
    end
    
end
%% load or calculate grid
if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefile,'mnu4Sq','sin2T4','chi2','chi2_ref','FitResults_Null')
    fprintf('load grid from file %s \n',savefile)
else
    if range<=40
        FSDFlag = 'Sibille0p5eV';
    elseif range>40
        FSDFlag = 'SibilleFull';
    end
    
    RunAnaArg = {'RunNr',1,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag',FSDFlag,...
        'ELossFlag',ELossFlag,...
        'AngularTFFlag',AngularTFFlag,...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'TwinBias_Q',18573.73,...
        'SysBudget',SysBudget,...
        'pullFlag',pullFlag,...
        'FakeInitFile',FakeInitFile,...
        'PixList',1:136};
    
    T = RunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    T.chi2 = chi2;
    if ~strcmp(chi2,'chi2Stat')
        AllEffects = struct(...
            'RF_EL','OFF',...   % Response Function(RF) EnergyLoss
            'RF_BF','ON',...   % RF B-Fields
            'RF_RX','ON',...   % Column Density, inel cross ection
            'FSD','ON',...
            'TASR','OFF',...
            'TCoff_RAD','OFF',...
            'TCoff_OTHER','ON',...
            'DOPoff','OFF',...
            'Stack','OFF',...
            'FPDeff','OFF',...
            'LongPlasma','OFF');
        T.ComputeCM('SysEffect',AllEffects,'BkgCM','OFF');
    end
    %% null hypothesis : no steriles
    T.Fit;
    FitResults_Null = T.FitResult;
    
    %% define msq4 - sin2t4 grid
    sin2T4      = logspace(log10(4e-04),log10(0.5),nGridSteps); %linspace(0.001,0.5,nGridSteps)
    mnu4Sq      = logspace(-1,log10((range+5)^2),nGridSteps)';
    mnu4Sq      = repmat(mnu4Sq,1,nGridSteps);
    sin2T4      = repmat(sin2T4,nGridSteps,1);
    %% make copy of model for parallel computing
    D = copy(repmat(T,nGridSteps.*nGridSteps,1));
    
    %% scan over msq4-sin2t4 grid
    D              = reshape(D,numel(D),1);
    chi2Grid       = zeros(nGridSteps*nGridSteps,1);
    FitResultsGrid = cell(nGridSteps*nGridSteps,1);
    mnu4Sq_Grid    = reshape(mnu4Sq',nGridSteps*nGridSteps,1);
    sin2T4_Grid    = reshape(sin2T4',nGridSteps*nGridSteps,1);
    
    parfor i= 1:(nGridSteps*nGridSteps)
        D(i).SimulateRun;
        D(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
        D(i).Fit
        chi2Grid(i) = D(i).FitResult.chi2min;
        FitResultsGrid{i} = D(i).FitResult;
    end
    
    chi2       = reshape(chi2Grid,nGridSteps,nGridSteps);
    FitResults = reshape(FitResultsGrid,nGridSteps,nGridSteps);
    
    if min(min(chi2))<T.FitResult.chi2min
        chi2_ref = min(min(chi2));
    else
        chi2_ref = T.FitResult.chi2min;
    end
    %% save
    save(savefile,'chi2_ref','RunAnaArg',...%'FitResults_ref'
        'chi2','mnu4Sq','sin2T4','FitResults','FitResults_Null');
    fprintf('save file to %s \n',savefile)
end
end


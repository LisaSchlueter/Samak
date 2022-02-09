% Choose random pixellist with half of KNM1 pixels
savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
nPixList = 1000;
FitResults = cell(nPixList,1);
PixLists   = cell(nPixList,1);
chi2 = 'chi2Stat';
NP = 1.064;

savename = [savedir,sprintf('RandomHalfPixList_Unblinded_%s_NP%2g_%.0ffits.mat',chi2,NP,nPixList)];
%%
if exist(savename,'file')
    load(savename);
else
    MakeDir(savedir);
    
    RunAnaArg = {'RunList','KNM1',...
        'chi2','chi2Stat',...
        'DataType','Real',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',NP,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF'};
    
    M = MultiRunAnalysis(RunAnaArg{:});
    M.Fit;
    PixList_def = M.PixList; % regular golden pixel list
    nPix = numel(M.PixList);
    nPixHalf     = ceil(nPix/2);           % take only half of all runs
    
    D = copy(repmat(M,nPixList,1));
    D = reshape(D,numel(D),1);
    %%
    
    parfor i=1:nPixList
        RandIndex    = randperm(nPix);         % randomly permute runlist indices
        PixLists{i}   = PixList_def(RandIndex(1:nPixHalf));
        
        D(i).PixList =  PixLists{i};
        D(i).StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
        
        D(i).SimulateStackRuns;
        D(i).exclDataStart = D(i).GetexclDataStart(40);
        D(i).Fit;
        
        
        FitResults{i} = D(i).FitResult;
        PixLists{i}   = D(i).RunList; % wrong....ups, not so important though
      
    end
    
    save(savename,'FitResults','PixLists','M','RunAnaArg');
    fprintf('save to %s \n',savename);
end

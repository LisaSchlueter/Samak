% Uniform fit on KNM2 data
% random half of all golden runs
% March 2022, Lisa

nRunList = 2000;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;               % fit range in eV below endpoint
FSDFlag = 'KNM2_0p1eV';
BKG_PtSlope = 3*1e-06;
%savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savename = sprintf('%sknm2_RunListRandHalf_%s_%s_%.0feV_%.0ffits_%s_BkgPT%.2gmucpsPers.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nRunList,FSDFlag,BKG_PtSlope*1e6);


if exist(savename,'file')
    load(savename,'RunList','FitResult','RunAnaArg')
else
    MakeDir(savedir);
    
    SigmaSq =  0.0124+0.0025;
    RunList = cell(nRunList,1);
    FitResult = cell(nRunList,1);
    RunAnaArg = {'RunList','KNM2_RandHalf',...  % define run number -> see GetRunList
        'fixPar',freePar,...         % free Parameter !!
        'DataType',DataType,...              % Real, Twin or Fake
        'FSDFlag',FSDFlag,...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1.112,...
        'MosCorrFlag','OFF',...
        'TwinBias_Q',18573.7,...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'BKG_PtSlope',BKG_PtSlope};
    
    M = MultiRunAnalysis(RunAnaArg{:});
    M.exclDataStart = M.GetexclDataStart(range);
    
    D = copy(repmat(M,nRunList,1));
    D = reshape(D,numel(D),1);
    FitResults = cell(nRunList,1);
    RunLists   = cell(nRunList,1);
    
    parfor i=1:nRunList
        %% build object of MultiRunAnalysis class
        D(i).RunList = D(i).GetRunList;
        D(i).nRuns         = length(D(i).RunList);
        D(i).StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
        D(i).SimulateStackRuns;
        
        D(i).exclDataStart = D(i).GetexclDataStart(40);
        D(i).Fit;

        FitResults{i} = D(i).FitResult;
        RunLists{i} = D(i).RunList;
    end
    
    save(savename,'FitResults','RunLists','M','RunAnaArg');
    fprintf('save to %s \n',savename);
end

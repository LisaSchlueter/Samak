% Fit KNM2 single pixels 

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
chi2 = 'chi2Stat';
NP = 1.064;

savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
savename = sprintf('%sknm1_SinglePixFit_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,chi2,NP);
%%
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    RunAnaArg = {'RunList','KNM1',...
        'RingMerge','None',...
        'chi2',chi2,...
        'DataType','Real',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',NP,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF'};
    
    M = MultiRunAnalysis(RunAnaArg{:});          % init model, read data
    M.exclDataStart = M.GetexclDataStart(range); % set fit range
    M.Fit;
    FitResult_Uniform = M.FitResult;
    PixList = M.PixList;
    
    FitResult = cell(148,1);
    mNuSq  = NaN.*zeros(148,1);
    mNuSqErr  = NaN.*zeros(148,1);
    
    
    parfor i=1:numel(PixList)
        
        A = MultiRunAnalysis(RunAnaArg{:},'PixList',PixList(i));          % init model, read data
        A.exclDataStart = A.GetexclDataStart(40); % set fit range
        A.Fit;
        FitResult{i} = A.FitResult;
        mNuSq(i) = A.FitResult.par(1);
        mNuSqErr(i) = 0.5.*(A.FitResult.errPos(1)-A.FitResult.errNeg(1));
    end
    
    save(savename,'M','PixList','FitResult_Uniform',...
        'FitResult','mNuSq','mNuSqErr');
end





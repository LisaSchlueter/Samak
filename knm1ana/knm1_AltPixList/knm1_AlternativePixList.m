% Fit KNM2 alternative runlist

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
chi2 = 'chi2Stat';
NP = 1.064;

AltPixLists = {'AziHalfNS','AziHalfEW'};%,'Half'};
RecomputeFlag = 'OFF';
for i=1:numel(AltPixLists)
    AltPixList =AltPixLists{i};  %'Half';%'AziHalfEW';  % defines alternative pixel list
    % label
    savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
    savename = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
        savedir,AltPixList,DataType,strrep(freePar,' ',''),range,chi2,NP);
    
    if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
        load(savename);
    else
        RunAnaArg = {'RunList','KNM1',...
            'RingMerge',AltPixList,...
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
        
        R = RingAnalysis('RunAnaObj',M,'RingList',M.RingList);
        R.FitRings;
        
        FitResult = R.FitResult(1);
        PixList = arrayfun(@(x) x.PixList,R.MultiObj,'UniformOutput',false);
        RunList = M.RunList;
        Q_i = M.ModelObj.Q_i;
        
        mNuSq    = FitResult.par(:,1);
        mNuSqErr = FitResult.err(:,1);
        E0       = FitResult.par(:,2)+Q_i;
        E0Err    = FitResult.err(:,2);
        
        MakeDir(savedir);
        save(savename,'FitResult','PixList','RunList','Q_i','mNuSq','E0','mNuSqErr','E0Err','RunAnaArg');
    end
end

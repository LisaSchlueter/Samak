% single run fit
% stat + NP only
% nu-mass fix
range   = 40;
freePar = 'E0 Bkg Norm';
chi2    = 'chi2Stat';
SysBudget = 36;

savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/'];
savename = sprintf('%sknm2ub1_FitRunList_%.0feV_%s_%s.mat',...
    savedir,range,strrep(freePar,' ',''),chi2);

if ~strcmp(chi2,'chi2Stat')
    savename = strrep(savename,'.mat',sprintf('_SysBudget%s.mat',SysBudget));
end

if exist(savename,'file')
    load(savename,'FitResult','RunAnaArg','A');
else
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2','chi2Stat',...
        'DataType','Real',...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag','StackPixel',...
        'RingMerge','Full',...
        'chi2',chi2,...
        'pullFlag',99,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1.112};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);

    FitResult = A.FitRunList;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A')
end

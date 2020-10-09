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
        'NonPoissonScaleFactor',1};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);

    FitResult = A.FitRunList;
    MakeDir(savedir);
    save(savename,'FitResult','RunAnaArg','A')
end

%% plot
saveplot = 'pdf';
FitResults= A.PlotFitRunList('Parameter','B','YLim',[170 280],'saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','E0','YLim',[-0.7,1],'DisplayStyle','Rel','saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','N','YLim',[0.9,1.1],'saveplot',saveplot,'HideGaps','OFF');
A.PlotFitRunList('Parameter','pVal','YLim',[-0.2,1.2],'saveplot',saveplot,'HideGaps','OFF');

%%
E0 = A.SingleRun_FitResults.chi2Stat.E0;
B = A.SingleRun_FitResults.chi2Stat.B;
N = A.SingleRun_FitResults.chi2Stat.N;
p = A.SingleRun_FitResults.chi2Stat.pValue;
fprintf('--------------------------------------\n')
fprintf('Run list: %s (%.0f runs) \n',A.RunData.RunName,A.nRuns)
fprintf('<E0> = %.3f eV , std = %.3f eV \n',mean(E0),std(E0));
fprintf('<B>  = %.3f mcps , std = %.3f mcps \n',1e3.*mean(B),1e3.*std(B));
fprintf('<N>  = %.3f       , std = %.3f  \n',1+mean(N),std(N));
fprintf('<pval> = %.2f     , std = %.2f   \n ',mean(p),std(p));
fprintf('--------------------------------------\n')

close all

%% single period
Period = 2;

switch Period
    case 1
        nRuns = 171;
        Start = 1;
        Stop = nRuns;
    case 2
        nRuns = 97;
        Start = 172;
        Stop = Start+nRuns-1;
    case 3
        nRuns = 93;
        Start = 268;
        Stop = 360;
end
fprintf('--------------------------------------\n')
fprintf('Run list Period: %.0f (%.0f runs) \n',Period,nRuns)
fprintf('<E0> = %.3f eV , std = %.3f eV \n',mean(E0(Start:Stop)),std(E0(Start:Stop)));
fprintf('<B>  = %.1f mcps , std = %.1f mcps \n',1e3.*mean(B(Start:Stop)),1e3.*std(B(Start:Stop)));
fprintf('<N>  = %.3f       , std = %.3f  \n',1+mean(N(Start:Stop)),std(N(Start:Stop)));
fprintf('<pval> = %.2f     , std = %.2f   \n ',mean(p(Start:Stop)),std(p(Start:Stop)));
fprintf('--------------------------------------\n')

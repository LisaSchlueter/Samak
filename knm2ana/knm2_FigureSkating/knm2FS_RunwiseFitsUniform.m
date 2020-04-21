% KNM2 Figure skating twins

range = 40;

RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...              
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'ROIFlag','Default',...
    'SynchrotronFlag','ON',...
    'AngularTFFlag','ON'};

%% build object of MultiRunAnalysis class
D = MultiRunAnalysis(RunAnaArg{:});
D.exclDataStart = D.GetexclDataStart(range);

%%
saveplot = 'pdf';
FitResults= D.PlotFitRunList('Parameter','B','YLim',[170 280],'saveplot',saveplot,'HideGaps','OFF');
% D.PlotFitRunList('Parameter','E0','YLim',[-0.7,1],'DisplayStyle','Rel','saveplot',saveplot,'HideGaps','OFF');
% D.PlotFitRunList('Parameter','N','YLim',[0.9,1.1],'saveplot',saveplot,'HideGaps','OFF');
% D.PlotFitRunList('Parameter','pVal','YLim',[-0.2,1.2],'saveplot',saveplot,'HideGaps','OFF');

%%
E0 = D.SingleRun_FitResults.chi2Stat.E0;
B = D.SingleRun_FitResults.chi2Stat.B;
N = D.SingleRun_FitResults.chi2Stat.N;
p = D.SingleRun_FitResults.chi2Stat.pValue;
fprintf('--------------------------------------\n')
fprintf('Run list: %s (%.0f runs) \n',D.RunData.RunName,D.nRuns)
fprintf('<E0> = %.3f eV , std = %.3f eV \n',mean(E0),std(E0));
fprintf('<B>  = %.3f mcps , std = %.3f mcps \n',1e3.*mean(B),1e3.*std(B));
fprintf('<N>  = %.3f       , std = %.3f  \n',1+mean(N),std(N));
fprintf('<pval> = %.2f     , std = %.2f   \n ',mean(p),std(p));
fprintf('--------------------------------------\n')

close all
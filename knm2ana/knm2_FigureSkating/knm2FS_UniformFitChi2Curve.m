% KNM2 Figure skating twins
range = 40;
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
RecomputeFlag = 'OFF';
chi2 = 'chi2Stat';%CMShape';
%% load or calc

savedir = [getenv('SamakPath'),'knm2ana/knm2_FigureSkating/results/'];
savename = sprintf('%sknm2FS_UniformFitChi2profile_%.0feV_%s.mat',savedir,range,chi2);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
        'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2',chi2,...              % statistics only
        'NonPoissonScaleFactor',1,...
        'TwinBias_Q',E0,...
        'ROIFlag','14keV',...
        'SysBudget',33};
    
    %% build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);

    if ~strcmp(chi2,'chi2Stat')
        A.NonPoissonScaleFactor = 1.112;
    end
    %% broaden RF +  broadening/shift of FSD
    TimeSec = zeros(3,1);
    TimeSec(1) = sum(A.SingleRunData.TimeSec(1:171));
    TimeSec(2) = sum(A.SingleRunData.TimeSec(172:268));
    TimeSec(3) = sum(A.SingleRunData.TimeSec(269:361));
    MultiWeights = TimeSec./sum(TimeSec);
    MultiPos = [E0(1),E0(end-120),E0(end)]';
    MultiPosRel = MultiPos-wmean(MultiPos,MultiWeights);
    Sigma = std(E0).^2;
    A.ModelObj.LoadFSD('MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
        'SanityPlot','ON','Sigma',Sigma);
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    MACE_Sigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = MACE_Sigma;
    A.ModelObj.InitializeRF;
    %%
    A.Fit;
    FitResult_imp  = A.FitResult;
   
    ScanResults = A.GetAsymFitError('Parameter','mNu',...
                                     'Mode','Uniform',...
                                     'ParScanMax',0.6,...
                                     'nFitMax',20);
    save(savename,'FitResult_imp','E0','MACE_Sigma','A','ScanResults');
end
%%
close all
outCM = A.PlotChi2Curve('FitResult',FitResult_imp,'ScanResult',ScanResults,...
    'Parameter','mNu','HoldOn','OFF');

plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(gcf,plotname);
%% if both stat and stat + syst run this (execute upper part twice: 1. chi2CMShape. 2. chi2Stat)
out = A.PlotChi2Curve('FitResult',FitResult_imp,'ScanResult',ScanResults,...
    'Parameter','mNu','HoldOn','ON');
out{2}.delete
out{3}.delete
out{4}.delete
outCM{2}.delete
outCM{3}.delete
outCM{4}.delete
leg = legend([out{1},outCM{1}],'Stat. only','Stat. and Syst');
leg.EdgeColor = rgb('Silver');
xlim([-0.4 0.4])
ylim([0 2])
 plotname = strrep(strrep(strrep(savename,'results','plots'),'.mat','.pdf'),chi2,'StatSyst');
 export_fig(gcf,plotname);
% KNM2 Figure skating twins
% Multi-Ring fit
% chi2 profile
% stat and stat + syst
RunList = 'KNM2_Prompt';
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
range = 40;
chi2 = 'chi2Stat';
pullFlag = 4;
freePar = 'mNu E0 Norm Bkg';
DataType = 'Twin';
RingMerge = 'Full';

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_MultiRingFitChi2profile_%s_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s.mat',...
    DataType, RunList,chi2,strrep(freePar,' ',''),pullFlag,range,RingMerge)];

if exist(savename,'file')
    load(savename)
else
    % read data and set up model
    RunArg = {'RunList',RunList,...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'SysBudget',33,...
        'AnaFlag','Ring',...
        'RingMerge',RingMerge,...
        'chi2','chi2Stat',...
        'pullFlag',pullFlag,...
        'TwinBias_Q',E0,...
        'ROIFlag','14keV',...
        'MosCorrFlag','OFF',...
        'NonPoissonScaleFactor',1};
    
    MR = MultiRunAnalysis(RunArg{:});
    MR.exclDataStart = MR.GetexclDataStart(range);
    MR.ModelObj.RFBinStep = 0.002;
    MR.ModelObj.InitializeRF;
    
    MR.Fit('SaveFit','OFF');
    FitResults_imp2 = MR.FitResult;
    
    if ~strcmp(chi2,'chi2Stat')
        MR.NonPoissonScaleFactor = 1.112;
        MR.SetNPfactor; % set correct dimensions
        MR.chi2 = chi2;
        MR.ComputeCM;
    end
    
    %% broaden of RF + broadening/shift of FSD
    TimeSec = zeros(3,1);
    TimeSec(1) = sum(MR.SingleRunData.TimeSec(1:171));
    TimeSec(2) = sum(MR.SingleRunData.TimeSec(172:268));
    TimeSec(3) = sum(MR.SingleRunData.TimeSec(269:361));
    MultiWeights = TimeSec./sum(TimeSec);
    MultiPos = [E0(1),E0(end-120),E0(end)]';
    MultiPosRel = repmat(MultiPos-wmean(MultiPos,MultiWeights),1,MR.nRings);
    Sigma = repmat(std(E0),3,MR.nRings);
    FSDArg = {'MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
        'SanityPlot','ON','Sigma',Sigma};
    MR.ModelObj.LoadFSD(FSDArg{:});
    
    %% fit and scan chi2 profile
    MR.Fit('SaveFit','OFF');
    FitResults_imp = MR.FitResult;
    
    ScanResults = MR.GetAsymFitError('Parameter','mNu',...
        'Mode','Uniform',...
        'ParScanMax',0.6,...
        'nFitMax',20);
    save(savename,'FitResults_imp','E0','MR','ScanResults','FSDArg');
end

%%
outCM = MR.PlotChi2Curve('FitResult',FitResults_imp,'ScanResult',ScanResults,...
    'Parameter','mNu','HoldOn','OFF');
%% if both stat and stat + syst run this (execute upper part twice: 1. chi2CMShape. 2. chi2Stat)
out = MR.PlotChi2Curve('FitResult',FitResults_imp,'ScanResult',ScanResults,...
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
fprintf('save plot to %s \n',plotname)
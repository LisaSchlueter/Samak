% KNM2 twin multi ring fit
% March 2020, Lisa

RunList = 'KNM2_Prompt';
E0 = knm2FS_GetE0Twins('SanityPlot','OFF');
range = 40;
chi2 = 'chi2CMShape';
pullFlag = 4;
freePar = 'mNu E0 Norm Bkg';
DataType = 'Twin';

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
    'RingMerge','Full',...
    'chi2','chi2Stat',...
    'pullFlag',pullFlag,...
    'TwinBias_Q',E0,...
    'ROIFlag','14keV',...
    'MosCorrFlag','OFF',...
    'NonPoissonScaleFactor',1};

MR = MultiRunAnalysis(RunArg{:});
%MR.ModelObj.recomputeRF ='ON';
%MR.ModelObj.InitializeRF;
%MR.ModelObj.recomputeRF ='OFF';
MR.exclDataStart = MR.GetexclDataStart(range);
if ~strcmp(chi2,'chi2Stat')
    MR.NonPoissonScaleFactor = 1.112;
    MR.chi2 = chi2;
    MR.ComputeCM;
end

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_MultiRingFit_%s_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s_01binning.mat',...
    DataType, RunList,chi2,strrep(freePar,' ',''),pullFlag,range,MR.RingMerge)];

if exist(savename,'file')
    load(savename)
    MR.FitResult = FitResults;
else
    MR.Fit('SaveFit','ON');
    FitResults = MR.FitResult;
    save(savename,'FitResults','RunArg');
    
    %% fit with broadening of RF + broadening/shift of FSD
    TimeSec = zeros(3,1);
    TimeSec(1) = sum(MR.SingleRunData.TimeSec(1:171));
    TimeSec(2) = sum(MR.SingleRunData.TimeSec(172:268));
    TimeSec(3) = sum(MR.SingleRunData.TimeSec(269:361));
    MultiWeights = TimeSec./sum(TimeSec);
    MultiPos = [E0(1),E0(end-120),E0(end)]';
    MultiPosRel = repmat(MultiPos-wmean(MultiPos,MultiWeights),1,MR.nRings);
    MR.ModelObj.recomputeRF='ON';
    MR.ModelObj.InitializeRF;
    Sigma = repmat(std(E0).^2,3,MR.nRings);
    MR.ModelObj.LoadFSD('MultiPos',MultiPosRel,'MultiWeight',MultiWeights,...
        'SanityPlot','ON','Sigma',Sigma);
    MR.Fit('SaveFit','OFF');
    FitResults_imp = MR.FitResult;
    save(savename,'FitResults_imp','-append');
end

% MR.PlotFitMultiRing('PlotPar','qU','linFit','ON','savePlot','ON','Blind','OFF');
%
% plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',MR.DataSet)];
% freePar = ConvertFixPar('freePar',MR.fixPar,'nPar',MR.nPar,'nPixels',numel(MR.RunData.MACE_Ba_T),'Mode','Reverse');
% plotname = [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_freePar%s_%s.pdf',...
%     MR.RingMerge,MR.RunData.RunName,MR.chi2,freePar,'qU')];
% PlotRingWiseFitPar(MR,'SaveAs',plotname); %FPD viewer
% 

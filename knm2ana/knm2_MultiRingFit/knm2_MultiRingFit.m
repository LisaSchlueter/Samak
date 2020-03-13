% Script to develop and test multi ring fit
% based on KNM2 runs
% October 2019, Lisa

RunList = 'KNM2_RW3';%'KNM2_afterFix';%'KNM2_RW3';%

%exclDataStart= 11;
range = 40;
chi2 = 'chi2Stat';%CMShape';
pullFlag = 4;
fixPar = 'E0 Norm Bkg qU';

% read data and set up model
RunArg = {'RunList',RunList,...
    'chi2','chi2Stat','DataType','Twin',...
    'fixPar',fixPar,...
    'RadiativeFlag','ON',...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AnaFlag','Ring',...
    'RingMerge','Full',...
    'chi2',chi2,...
    'pullFlag',pullFlag,...
    'TwinBias_Q',18573.70,...
    'ROIFlag','Default',...
    'MosCorrFlag','OFF'};

MR = MultiRunAnalysis(RunArg{:});
MR.exclDataStart = MR.GetexclDataStart(range);

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_MultiRingFit_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s.mat',...
    RunList,chi2,strrep(fixPar,' ',''),pullFlag,range,MR.RingMerge)];

if exist(savename,'file')
    load(savename)
    MR.FitResult = FitResults;
else
    MR.Fit('SaveFit','ON');
    FitResults = MR.FitResult;
    save(savename,'FitResults','RunArg');
end
%%
MR.PlotFitMultiRing('PlotPar','qU','linFit','ON','savePlot','ON','Blind','OFF');

plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',MR.DataSet)];
freePar = ConvertFixPar('freePar',MR.fixPar,'nPar',MR.nPar,'nPixels',numel(MR.RunData.MACE_Ba_T),'Mode','Reverse');
plotname = [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_freePar%s_%s.pdf',...
    MR.RingMerge,MR.RunData.RunName,MR.chi2,freePar,'qU')];
PlotRingWiseFitPar(MR,'SaveAs',plotname); %FPD viewer


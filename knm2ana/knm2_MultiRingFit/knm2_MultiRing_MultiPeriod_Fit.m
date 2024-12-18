% Script to develop and test multi ring + multi period fit
% based on KNM2 runs
% March 2020, Lisa
tic;

%% settings
RunList  = 'KNM2_Prompt'; % all data -> multi-period
range    = 40;
chi2     = 'chi2Stat';%CMShape';
pullFlag = 4;
fixPar   = 'E0 Norm Bkg';

%% twins with plasma drift
E0OffseteV  = [0,0.1,-0.1]';          % per RW-perid (the same for all rings )
DriftPerDay = [6*1e-03, 0, 6*1e-03]'; % per RW-perid (the same for all rings )

[E0,RectWidth,MultiWeights] = ConvertPlasmaDriftDay2E0Run('DriftPerDay',DriftPerDay,...
    'E0OffseteV',E0OffseteV,...
    'E0ref',18573.70,...
    'SanityPlot','OFF');

MultiPos   =  repmat(E0OffseteV,[1,4]); % size: nPeriod x nPseudo-Rings --> here eg. 3x4
MultiSigma =  repmat(RectWidth,[1,4]);
MultiWeights = [1/3 1/3 1/3];
%% read data and set up model: twins with plasma drift and steps
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
    'ROIFlag','14keV',...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',E0};

MR = MultiRunAnalysis(RunArg{:});
MR.exclDataStart = MR.GetexclDataStart(range);

%% Modify FSD
MR.ModelObj.LoadFSD('MultiPos',MultiPos,...
                    'Sigma',MultiSigma,...
                    'MultiWeights',MultiWeights,...
                    'SanityPlot','ON');
                
MR.ModelObj.ComputeTBDDS;
MR.ModelObj.ComputeTBDIS;

%% label 
savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_MultiRingFit_MultiPeriod_%s_%s_%s_pull%.0f_%.0feVrange_RingMerge%s.mat',...
             RunList,chi2,strrep(fixPar,' ',''),pullFlag,range,MR.RingMerge)];
%%
if exist(savename,'file')
    load(savename)
    MR.FitResult = FitResults;
else
    MR.Fit('SaveFit','ON');
    FitResults = MR.FitResult;
    save(savename,'FitResults','RunArg');
end
%% display fit result (optional)
MR.PlotFitMultiRing('PlotPar','qU','linFit','ON','savePlot','ON','Blind','OFF');

plotdir = [getenv('SamakPath'),sprintf('tritium-data/plots/%s/MultiRingFit/',MR.DataSet)];
freePar = ConvertFixPar('freePar',MR.fixPar,'nPar',MR.nPar,'nPixels',numel(MR.RunData.MACE_Ba_T),'Mode','Reverse');
plotname = [plotdir,sprintf('FPDViewer_MultiRing%s_%s_%s_freePar%s_%s.pdf',...
    MR.RingMerge,MR.RunData.RunName,MR.chi2,freePar,'qU')];
PlotRingWiseFitPar(MR,'SaveAs',plotname); %FPD viewer
toc;

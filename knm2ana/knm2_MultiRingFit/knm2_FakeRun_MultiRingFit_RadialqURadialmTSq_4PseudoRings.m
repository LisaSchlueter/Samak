% Script to develop and test multi ring fit
% based on KNM2 fake MC run
% October 2019, Lisa

% method: 
% test qU-Offset:
% Data: MC data with radial dependent qU-values. 2 rings with offset (defined in InitFile)
% Model: Test different cases:
%        1) Uniform fit with averaged qU's. Expect a neutrino mass shift
%        2) Multi ring fit with correct ringwise qU-values
%        3) Multi ring fit with uniform qU-values. Try to fit qU-Offset
%% init file
InitFile = @ref_FakeRun_KNM2_CD84_50days_RadialqU; %100 column density, 50 days
Display = 'ON';

%% set up model, compute fake run if necessary
RunAnaArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',11,... % 11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'RingMerge','Full',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1,...
    'FakeInitFile',InitFile,...
    'AnaFlag','Ring',...
    'fixPar','mNu E0 Bkg Norm qU mTSq',...
    'PullFlag',4,'fitter','minuit'};
M = RunAnalysis(RunAnaArg{:});

mTSqBias =[0, 0.1,0.2,0.3];

% label
savedirPlot  = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/plots/'];
savenamePlot = [savedirPlot,sprintf('knm2_FakeRun_MultiRingFit_RadialqUfit_4Rings_mTSqfree_Sigma2_%.1f-%.1feV2.pdf',min(mTSqBias),max(mTSqBias))];
savedir = strrep(savedirPlot,'plots','results');
MakeDir(savedir); MakeDir(savedirPlot);
savenameResult = [savedir,sprintf('knm2_FakeRun_MultiRingFit_RadialqUfit_4Rings_mTSqfree_Sigma2_%.1f-%.1feV2.mat',min(mTSqBias),max(mTSqBias))];

if exist(savenameResult,'file')
    load(savenameResult);
    M.FitResult = FitResult;
else
%%  multi-ring fit, take qU-values from first ring
M.ModelObj.ComputeTBDDS('mTSq_bias',mTSqBias);
M.ModelObj.ComputeTBDIS;
TBDIS = M.ModelObj.TBDIS;

qUfirstring = repmat(M.RunData.qU(:,1),[1,4]);
M.SimulateRun('qU',qUfirstring);

M.RunData.TBDIS = TBDIS;
M.RunData.TBDISE = sqrt(TBDIS);
M.Fit;
if strcmp(Display,'ON')
    [plotHandle, cbHandle]  = PlotRingWiseFitPar(M,'PlotPar','qU');
    export_fig(plotHandle,savenamePlot);
end

FitResult = M.FitResult;
save(savenameResult,'FitResult','RunAnaArg','mTSqBias','-mat');
end

%%
M.PlotFitMultiRing('PlotPar','mTSq','linFit','OFF','savePlot','ON')
